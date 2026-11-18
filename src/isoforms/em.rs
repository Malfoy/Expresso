use super::{InferArgs, Model, counts, load};
use crate::{abundance, export, index, input, output};
use anyhow::{Context, Result, ensure};
use rayon::prelude::*;
use serde::Serialize;
use serde_json::json;
use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fs,
    path::{Path, PathBuf},
};

struct Group {
    members: Vec<usize>,
    length: u64,
    profile: Vec<(usize, u64)>,
}
struct Row {
    feature: usize,
    terms: Vec<(usize, f64)>,
}
struct Component {
    groups: Vec<usize>,
    rows: Vec<Row>,
    rank: Option<usize>,
}
struct Prepared {
    groups: Vec<Group>,
    components: Vec<Component>,
    transcript_group: Vec<usize>,
    group_component: Vec<usize>,
}

// Sparse column elimination on columns scaled to unit maximum; rank is numerical.
fn rank(columns: impl Iterator<Item = Vec<(usize, f64)>>) -> usize {
    let mut basis: BTreeMap<usize, BTreeMap<usize, f64>> = BTreeMap::new();
    for column in columns {
        let scale = column.iter().map(|(_, v)| v.abs()).fold(0.0, f64::max);
        if scale == 0.0 {
            continue;
        }
        let mut column: BTreeMap<_, _> = column
            .into_iter()
            .map(|(row, v)| (row, v / scale))
            .collect();
        loop {
            column.retain(|_, v| v.abs() > 1e-10);
            let Some((&pivot, &coefficient)) = column.first_key_value() else {
                break;
            };
            if let Some(old) = basis.get(&pivot) {
                column.remove(&pivot);
                for (&row, &value) in
                    old.range((std::ops::Bound::Excluded(pivot), std::ops::Bound::Unbounded))
                {
                    *column.entry(row).or_default() -= coefficient * value;
                }
            } else {
                for value in column.values_mut() {
                    *value /= coefficient;
                }
                basis.insert(pivot, column);
                break;
            }
        }
    }
    basis.len()
}
fn root(parents: &mut [usize], id: usize) -> usize {
    let mut top = id;
    while parents[top] != top {
        top = parents[top];
    }
    let mut current = id;
    while parents[current] != current {
        let next = parents[current];
        parents[current] = top;
        current = next;
    }
    top
}
impl Prepared {
    fn new(model: &Model, rank_limit: usize) -> Result<Self> {
        let mut profiles = vec![Vec::new(); model.transcripts.len()];
        for (id, feature) in model.features.iter().enumerate() {
            for occurrence in &feature.transcripts {
                profiles[occurrence.member].push((
                    id,
                    occurrence
                        .multiplicity
                        .checked_mul(feature.distinct_kmers)
                        .context("exposure overflow")?,
                ));
            }
        }
        let mut groups: Vec<Group> = Vec::new();
        let mut intern: HashMap<Vec<(usize, u64)>, usize> = HashMap::new();
        let mut transcript_group = Vec::new();
        for (id, profile) in profiles.into_iter().enumerate() {
            let existing = if profile.is_empty() {
                None
            } else {
                intern.get(&profile).copied()
            };
            let group = if let Some(group) = existing {
                groups[group].members.push(id);
                group
            } else {
                let group = groups.len();
                if !profile.is_empty() {
                    intern.insert(profile.clone(), group);
                }
                groups.push(Group {
                    members: vec![id],
                    length: model.transcripts[id].observable_kmers,
                    profile,
                });
                group
            };
            transcript_group.push(group);
        }
        drop(intern);
        let mut parents: Vec<_> = (0..groups.len()).collect();
        let mut rows = Vec::new();
        for (feature_id, feature) in model.features.iter().enumerate() {
            let mut terms = BTreeMap::new();
            for m in &feature.transcripts {
                let group = transcript_group[m.member];
                let value = m.multiplicity * feature.distinct_kmers;
                if let Some(old) = terms.insert(group, value) {
                    ensure!(old == value, "inconsistent identical-profile group");
                }
            }
            let first = *terms.first_key_value().unwrap().0;
            for &id in terms.keys().skip(1) {
                let a = root(&mut parents, first);
                let b = root(&mut parents, id);
                if a != b {
                    parents[b] = a;
                }
            }
            rows.push((feature_id, terms));
        }
        let mut components: Vec<Component> = Vec::new();
        let mut component_ids = BTreeMap::new();
        let mut group_component = Vec::new();
        let mut local = Vec::new();
        for group in 0..groups.len() {
            let top = root(&mut parents, group);
            let id = *component_ids.entry(top).or_insert_with(|| {
                components.push(Component {
                    groups: Vec::new(),
                    rows: Vec::new(),
                    rank: None,
                });
                components.len() - 1
            });
            group_component.push(id);
            local.push(components[id].groups.len());
            components[id].groups.push(group);
        }
        for (feature, terms) in rows {
            let id = group_component[*terms.first_key_value().unwrap().0];
            components[id].rows.push(Row {
                feature,
                terms: terms
                    .into_iter()
                    .map(|(g, value)| (local[g], value as f64))
                    .collect(),
            });
        }
        for component in &mut components {
            if component.groups.len() <= rank_limit {
                component.rank = Some(rank(component.groups.iter().map(|&g| {
                    groups[g]
                        .profile
                        .iter()
                        .map(|&(r, v)| (r, v as f64))
                        .collect()
                })));
            }
        }
        Ok(Self {
            groups,
            components,
            transcript_group,
            group_component,
        })
    }
}

#[derive(Serialize)]
struct Diagnostic {
    iterations: usize,
    converged: bool,
    active_features: usize,
    observed_sum: f64,
    log_likelihood: f64,
    fixed_point_change: f64,
    relative_l1_residual: f64,
}
struct ComponentFit {
    rates: Vec<f64>,
    assigned: Vec<f64>,
    diagnostic: Diagnostic,
}
struct Fit {
    rates: Vec<f64>,
    assigned: Vec<f64>,
    diagnostics: Vec<Diagnostic>,
}
fn expectation(
    component: &Component,
    y: &[f64],
    rates: &[f64],
    assigned: &mut [f64],
) -> Result<f64> {
    assigned.fill(0.0);
    let mut likelihood = 0.0;
    for row in &component.rows {
        let count = y[row.feature];
        if count == 0.0 {
            continue;
        }
        let lambda: f64 = row.terms.iter().map(|&(t, a)| a * rates[t]).sum();
        ensure!(
            lambda.is_finite() && lambda > 0.0,
            "positive feature has zero/nonfinite predicted abundance; check input magnitude"
        );
        likelihood += count * lambda.ln();
        for &(t, a) in &row.terms {
            assigned[t] += count * (a * rates[t] / lambda);
        }
    }
    Ok(likelihood)
}
fn fit_component(
    component: &Component,
    groups: &[Group],
    y: &[f64],
    max_iterations: usize,
    tolerance: f64,
) -> Result<ComponentFit> {
    let lengths: Vec<_> = component
        .groups
        .iter()
        .map(|&g| groups[g].length as f64)
        .collect();
    let observed_sum: f64 = component.rows.iter().map(|r| y[r.feature]).sum();
    ensure!(
        observed_sum.is_finite(),
        "feature total exceeds floating-point range"
    );
    let exposure: f64 = lengths.iter().sum();
    let initial = if exposure > 0.0 {
        observed_sum / exposure
    } else {
        0.0
    };
    let mut rates = vec![initial; component.groups.len()];
    let mut assigned = vec![0.0; rates.len()];
    let mut iterations = 0;
    let mut previous_likelihood = f64::NEG_INFINITY;
    if observed_sum > 0.0 {
        for iteration in 1..=max_iterations {
            let likelihood = expectation(component, y, &rates, &mut assigned)?
                - rates.iter().zip(&lengths).map(|(r, l)| r * l).sum::<f64>();
            ensure!(likelihood.is_finite(), "nonfinite EM likelihood");
            let mut change = 0.0;
            let mass: f64 = rates.iter().sum();
            for (t, value) in rates.iter_mut().enumerate() {
                let next = assigned[t] / lengths[t];
                ensure!(next.is_finite() && next >= 0.0, "invalid EM update");
                change += (next - *value).abs();
                *value = next;
            }
            iterations = iteration;
            let relative = change / mass.max(f64::MIN_POSITIVE);
            if relative <= tolerance
                && (likelihood - previous_likelihood).abs() <= tolerance * observed_sum.max(1.0)
            {
                break;
            }
            previous_likelihood = likelihood;
        }
    }
    let likelihood = expectation(component, y, &rates, &mut assigned)?
        - rates.iter().zip(&lengths).map(|(r, l)| r * l).sum::<f64>();
    let fixed_point_change = if observed_sum > 0.0 {
        rates
            .iter()
            .enumerate()
            .map(|(t, &r)| (assigned[t] / lengths[t] - r).abs())
            .sum::<f64>()
            / rates.iter().sum::<f64>().max(f64::MIN_POSITIVE)
    } else {
        0.0
    };
    let residual: f64 = component
        .rows
        .iter()
        .map(|row| {
            let predicted: f64 = row.terms.iter().map(|&(t, a)| a * rates[t]).sum();
            (predicted - y[row.feature]).abs()
        })
        .sum();
    let active_features = component.rows.iter().filter(|r| y[r.feature] > 0.0).count();
    Ok(ComponentFit {
        rates,
        assigned,
        diagnostic: Diagnostic {
            iterations,
            converged: fixed_point_change <= tolerance,
            active_features,
            observed_sum,
            log_likelihood: likelihood,
            fixed_point_change,
            relative_l1_residual: if observed_sum > 0.0 {
                residual / observed_sum
            } else {
                0.0
            },
        },
    })
}
fn fit(prepared: &Prepared, y: &[f64], args: &InferArgs) -> Result<Fit> {
    let partials = prepared
        .components
        .par_iter()
        .map(|c| fit_component(c, &prepared.groups, y, args.max_iterations, args.tolerance))
        .collect::<Result<Vec<_>>>()?;
    let mut rates = vec![0.0; prepared.groups.len()];
    let mut assigned = vec![0.0; rates.len()];
    let mut diagnostics = Vec::new();
    for (component, result) in prepared.components.iter().zip(partials) {
        for (local, &g) in component.groups.iter().enumerate() {
            rates[g] = result.rates[local];
            assigned[g] = result.assigned[local];
        }
        diagnostics.push(result.diagnostic);
    }
    Ok(Fit {
        rates,
        assigned,
        diagnostics,
    })
}
fn status(prepared: &Prepared, fit: &Fit, group: usize) -> &'static str {
    let id = prepared.group_component[group];
    let component = &prepared.components[id];
    let diagnostic = &fit.diagnostics[id];
    if prepared.groups[group].length == 0 {
        "unobservable"
    } else if !diagnostic.converged {
        "non_converged"
    } else if diagnostic.observed_sum == 0.0 {
        "zero_evidence"
    } else if component.rank.is_none() {
        "rank_unchecked"
    } else if component.rank != Some(component.groups.len()) {
        "unresolved"
    } else if prepared.groups[group].members.len() > 1 {
        "group_only"
    } else {
        "identifiable"
    }
}
fn resolved(state: &str) -> bool {
    ["identifiable", "group_only", "zero_evidence"].contains(&state)
}

fn report(
    root: &Path,
    name: &str,
    model: &Model,
    prepared: &Prepared,
    fit: &Fit,
    compression: output::Compression,
    global: bool,
) -> Result<()> {
    let suffix = compression.suffix();
    let states: Vec<_> = (0..prepared.groups.len())
        .map(|g| status(prepared, fit, g))
        .collect();
    let total_rate: f64 = fit.rates.iter().sum();
    let ppm_available = states
        .iter()
        .all(|state| resolved(state) || *state == "unobservable");
    let mut gene_rates: BTreeMap<&str, (f64, bool)> = BTreeMap::new();
    for (g, group) in prepared.groups.iter().enumerate() {
        let genes: BTreeSet<_> = group
            .members
            .iter()
            .map(|&t| model.transcripts[t].gene.as_str())
            .collect();
        for gene in &genes {
            let entry = gene_rates.entry(gene).or_insert((0.0, true));
            if genes.len() == 1 {
                entry.0 += fit.rates[g];
            }
            entry.1 &= resolved(states[g]) && genes.len() == 1;
        }
    }
    let path = if global {
        root.join(format!("global{suffix}"))
    } else {
        root.join("datasets").join(format!("{name}{suffix}"))
    };
    output::csv(&path, compression, |csv| {
        csv.write_record([
            "transcript_id",
            "transcript_name",
            "gene_name",
            "group_id",
            "component_id",
            "observable_kmers",
            "abundance",
            "within_gene_fraction",
            "relative_abundance_ppm",
            "assigned_observations",
            "model_estimate",
            "status",
        ])?;
        for (t, transcript) in model.transcripts.iter().enumerate() {
            let g = prepared.transcript_group[t];
            let group = &prepared.groups[g];
            let state = states[g];
            let individual = group.members.len() == 1 && resolved(state);
            let rate = individual.then_some(fit.rates[g]);
            let gene = gene_rates[transcript.gene.as_str()];
            let fraction = rate
                .filter(|_| gene.1)
                .map(|r| if gene.0 > 0.0 { r / gene.0 } else { 0.0 });
            let ppm = rate.filter(|_| ppm_available).map(|r| {
                if total_rate > 0.0 {
                    1e6 * r / total_rate
                } else {
                    0.0
                }
            });
            let state = if group.members.len() > 1 && state == "zero_evidence" {
                "zero_evidence"
            } else if group.members.len() > 1 && resolved(state) {
                "group_only"
            } else {
                state
            };
            csv.serialize((
                t + 1,
                &transcript.name,
                &transcript.gene,
                g + 1,
                prepared.group_component[g] + 1,
                transcript.observable_kmers,
                rate,
                fraction,
                ppm,
                individual.then_some(fit.assigned[g]),
                (group.members.len() == 1 && group.length > 0).then_some(fit.rates[g]),
                state,
            ))?;
        }
        Ok(())
    })?;
    let group_path = root.join("groups").join(if global {
        PathBuf::from(format!("global{suffix}"))
    } else {
        PathBuf::from("datasets").join(format!("{name}{suffix}"))
    });
    output::csv(&group_path, compression, |csv| {
        csv.write_record([
            "group_id",
            "component_id",
            "gene_names",
            "transcript_names",
            "observable_kmers",
            "abundance",
            "relative_abundance_ppm",
            "assigned_observations",
            "model_estimate",
            "status",
        ])?;
        for (g, group) in prepared.groups.iter().enumerate() {
            let genes: BTreeSet<_> = group
                .members
                .iter()
                .map(|&t| model.transcripts[t].gene.as_str())
                .collect();
            let names = serde_json::to_string(
                &group
                    .members
                    .iter()
                    .map(|&t| &model.transcripts[t].name)
                    .collect::<Vec<_>>(),
            )?;
            let genes = serde_json::to_string(&genes)?;
            let rate = resolved(states[g]).then_some(fit.rates[g]);
            let ppm = rate.filter(|_| ppm_available).map(|r| {
                if total_rate > 0.0 {
                    1e6 * r / total_rate
                } else {
                    0.0
                }
            });
            csv.serialize((
                g + 1,
                prepared.group_component[g] + 1,
                genes,
                names,
                group.length,
                rate,
                ppm,
                fit.assigned[g],
                (group.length > 0).then_some(fit.rates[g]),
                states[g],
            ))?;
        }
        Ok(())
    })?;
    output::csv(
        &root.join("diagnostics").join(if global {
            PathBuf::from(format!("global{suffix}"))
        } else {
            PathBuf::from("datasets").join(format!("{name}{suffix}"))
        }),
        compression,
        |csv| {
            csv.write_record([
                "component_id",
                "groups",
                "features",
                "rank",
                "active_features",
                "observed_sum",
                "iterations",
                "converged",
                "log_likelihood",
                "fixed_point_change",
                "relative_l1_residual",
            ])?;
            for (id, (component, diagnostic)) in
                prepared.components.iter().zip(&fit.diagnostics).enumerate()
            {
                csv.serialize((
                    id + 1,
                    component.groups.len(),
                    component.rows.len(),
                    component.rank,
                    diagnostic.active_features,
                    diagnostic.observed_sum,
                    diagnostic.iterations,
                    diagnostic.converged,
                    diagnostic.log_likelihood,
                    diagnostic.fixed_point_change,
                    diagnostic.relative_l1_residual,
                ))?;
            }
            Ok(())
        },
    )
}
fn inputs(
    args: &InferArgs,
    model: &Model,
    hash: &[u8; 32],
) -> Result<(Vec<input::Dataset>, serde_json::Value)> {
    if let Some(fof) = &args.fof {
        return Ok((
            input::datasets(fof, abundance::Mode::Reads)?,
            json!({"complete":true,"source":"user-provided exact disjoint feature files; missing CSV rows mean zero"}),
        ));
    }
    let root = args.input.as_ref().unwrap().canonicalize()?;
    let manifest: serde_json::Value =
        serde_json::from_reader(output::reader(&root.join("manifest.json"))?)?;
    ensure!(
        manifest["format"] == "isoform-features"
            && manifest["version"] == 1
            && manifest["exact"] == true
            && manifest["model_sha256"] == crate::compact::hex(hash)
            && manifest["features"].as_u64() == Some(model.features.len() as u64),
        "expected exact feature results from this isoform model; ordinary exon/gene abundance vectors are insufficient"
    );
    let entries = manifest["datasets"]
        .as_array()
        .context("missing feature datasets")?;
    let mut datasets = Vec::new();
    let mut names = BTreeSet::new();
    for entry in entries {
        let name = export::text(entry, "name")?;
        ensure!(
            !name.is_empty()
                && name != "."
                && name != ".."
                && name
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b))
                && names.insert(name),
            "invalid or repeated dataset name"
        );
        datasets.push(input::Dataset {
            name: name.into(),
            path: export::source(&root, export::text(entry, "output")?)?,
            mode: abundance::Mode::Reads,
        });
    }
    ensure!(!datasets.is_empty(), "no feature datasets");
    Ok((
        datasets,
        json!({"coverage":manifest["coverage"],"failed_datasets":manifest["failed_datasets"]}),
    ))
}
pub fn infer(args: &InferArgs) -> Result<()> {
    ensure!(
        args.threads > 0
            && args.jobs > 0
            && args.max_iterations > 0
            && args.rank_limit > 0
            && args.tolerance.is_finite()
            && args.tolerance > 0.0
            && args.tolerance < 1.0,
        "invalid inference worker, iteration, rank or tolerance setting"
    );
    let (model, hash) = load(&args.index)?;
    let (datasets, coverage) = inputs(args, &model, &hash)?;
    let prepared = Prepared::new(&model, args.rank_limit)?;
    let stage = index::staging(&args.output)?;
    for dir in [
        "datasets",
        "groups",
        "groups/datasets",
        "diagnostics",
        "diagnostics/datasets",
    ] {
        fs::create_dir_all(stage.path().join(dir))?;
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    let mut totals = vec![0.0; prepared.groups.len()];
    let mut assigned = vec![0.0; prepared.groups.len()];
    let mut global_diagnostics: Vec<Diagnostic> = prepared
        .components
        .iter()
        .map(|_| Diagnostic {
            iterations: 0,
            converged: true,
            active_features: 0,
            observed_sum: 0.0,
            log_likelihood: 0.0,
            fixed_point_change: 0.0,
            relative_l1_residual: 0.0,
        })
        .collect();
    let mut entries = Vec::new();
    let jobs = args.jobs.min(args.threads).min(datasets.len());
    eprintln!(
        "Inferring {} datasets from {} disjoint features; {} transcripts, {} profile groups, {} independent components",
        datasets.len(),
        model.features.len(),
        model.transcripts.len(),
        prepared.groups.len(),
        prepared.components.len()
    );
    for batch in datasets.chunks(jobs) {
        let results = pool.install(|| {
            batch
                .par_iter()
                .map(|dataset| -> Result<Fit> {
                    let values = counts::read_counts(&dataset.path, model.features.len(), &hash)
                        .with_context(|| format!("dataset {}", dataset.name))?;
                    let fitted = fit(&prepared, &values, args)?;
                    report(
                        stage.path(),
                        &dataset.name,
                        &model,
                        &prepared,
                        &fitted,
                        args.compression,
                        false,
                    )?;
                    eprintln!(
                        "{}: isoform EM complete ({} non-converged components)",
                        dataset.name,
                        fitted.diagnostics.iter().filter(|d| !d.converged).count()
                    );
                    Ok(fitted)
                })
                .collect::<Result<Vec<_>>>()
        })?;
        for (dataset, fitted) in batch.iter().zip(results) {
            for (total, rate) in totals.iter_mut().zip(&fitted.rates) {
                *total += rate;
                ensure!(total.is_finite(), "global rate overflow");
            }
            for (total, value) in assigned.iter_mut().zip(&fitted.assigned) {
                *total += value;
            }
            for (global, diagnostic) in global_diagnostics.iter_mut().zip(&fitted.diagnostics) {
                global.iterations = global.iterations.max(diagnostic.iterations);
                global.converged &= diagnostic.converged;
                global.active_features = global.active_features.max(diagnostic.active_features);
                global.observed_sum += diagnostic.observed_sum;
                global.log_likelihood += diagnostic.log_likelihood;
                global.fixed_point_change =
                    global.fixed_point_change.max(diagnostic.fixed_point_change);
                global.relative_l1_residual +=
                    diagnostic.relative_l1_residual * diagnostic.observed_sum;
            }
            entries.push(json!({"name":dataset.name,"source":dataset.path,"output":format!("datasets/{}{}",dataset.name,args.compression.suffix()),
                "non_converged_components":fitted.diagnostics.iter().filter(|d| !d.converged).count()}));
        }
    }
    for diagnostic in &mut global_diagnostics {
        if diagnostic.observed_sum > 0.0 {
            diagnostic.relative_l1_residual /= diagnostic.observed_sum;
        }
    }
    let global = Fit {
        rates: totals,
        assigned,
        diagnostics: global_diagnostics,
    };
    report(
        stage.path(),
        "global",
        &model,
        &prepared,
        &global,
        args.compression,
        true,
    )?;
    let manifest = json!({"version":1,"format":"isoform-abundance","model_sha256":crate::compact::hex(&hash),"index":args.index.canonicalize()?,"datasets":entries,"source_coverage":coverage,
        "threads":args.threads,"jobs":jobs,"max_iterations":args.max_iterations,"tolerance":args.tolerance,"rank_limit":args.rank_limit,"rank_method":"sparse elimination on maximum-scaled exposure columns; tolerance 1e-10",
        "features":model.features.len(),"transcripts":model.transcripts.len(),"groups":prepared.groups.len(),"components":prepared.components.len(),
        "compression":args.compression,"global":format!("global{}",args.compression.suffix()),
        "semantics":"Poisson working model of disjoint exon/intersection counts; abundance is assigned observations / observable canonical-kmer positions. Identical profiles are grouped; unresolved/unobservable/non-converged individual abundances are blank. model_estimate is diagnostic, not an identified abundance. relative_abundance_ppm is among observable groups and blank when the total is unresolved; it is not fragment TPM. Global sums independently fitted rates and assigned observations, not per-dataset fractions or ppm. No confidence intervals are implied by this correlated-kmer working model."});
    fs::write(
        stage.path().join("manifest.json"),
        serde_json::to_vec_pretty(&manifest)?,
    )?;
    index::publish(stage, &args.output)?;
    eprintln!("Isoform abundance saved: {}", args.output.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn aggregate_exons_resolve_skipping_and_conserve_counts() {
        let groups = vec![
            Group {
                members: vec![0],
                length: 250,
                profile: vec![],
            },
            Group {
                members: vec![1],
                length: 200,
                profile: vec![],
            },
        ];
        let component = Component {
            groups: vec![0, 1],
            rank: Some(2),
            rows: vec![
                Row {
                    feature: 0,
                    terms: vec![(0, 100.0), (1, 100.0)],
                },
                Row {
                    feature: 1,
                    terms: vec![(0, 50.0)],
                },
                Row {
                    feature: 2,
                    terms: vec![(0, 100.0), (1, 100.0)],
                },
            ],
        };
        let result =
            fit_component(&component, &groups, &[1000.0, 300.0, 1000.0], 10000, 1e-10).unwrap();
        assert!(result.diagnostic.converged);
        assert!((result.rates[0] - 6.0).abs() < 1e-7 && (result.rates[1] - 4.0).abs() < 1e-7);
        assert!((result.assigned.iter().sum::<f64>() - 2300.0).abs() < 1e-9);
        assert!(result.diagnostic.relative_l1_residual < 1e-8);
        let zero = fit_component(&component, &groups, &[0.0; 3], 10, 1e-8).unwrap();
        assert_eq!(zero.rates, vec![0.0, 0.0]);
        let unfinished =
            fit_component(&component, &groups, &[1000.0, 300.0, 1000.0], 1, 1e-12).unwrap();
        assert!(!unfinished.diagnostic.converged);
    }
    #[test]
    fn rank_detects_dependencies_between_nonidentical_profiles() {
        // Distant alternative choices: AB + ab = Ab + aB, despite distinct columns.
        let columns = vec![
            vec![(0, 100.0), (2, 100.0)],
            vec![(0, 100.0), (3, 100.0)],
            vec![(1, 100.0), (2, 100.0)],
            vec![(1, 100.0), (3, 100.0)],
        ];
        assert_eq!(rank(columns.into_iter()), 3);
        assert_eq!(
            rank(vec![vec![(0, 1.0)], vec![(0, 1.0), (1, 1.0)]].into_iter()),
            2
        );
    }
}
