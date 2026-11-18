import {Distribution,parsePercentiles,formatStatistic,compareStatistics} from './statistics.js';
import {formatAbundance,normalize,foldChange} from './display.js';
import {Ranking,share,percentage,barplot} from './ranking.js';

// These views retain scalar summaries, a bounded quantile sample, and one visible
// profile page. Compressed abundance blocks remain managed by the single WASM worker.
export function createExplorer(context) {
  const {$}=context,PAGE=25;
  let overview,datasetResults,datasetPage=0,focus,focusPage=0;
  const catalog=()=>context.catalog(),integer=()=>context.normalization()==='none'&&catalog().encoding!=='efv';
  const realDatasets=(scope='all')=>catalog().datasets.flatMap((d,id)=>d.kind==='dataset'&&(scope==='all'||context.selectedDatasets().has(id))?[id]:[]);
  const targetIds=(scope='all')=>scope==='all'?Array.from({length:catalog().targets.length},(_,id)=>id):[...context.selectedTargets()].sort((a,b)=>a-b);
  const small=()=>catalog()&&catalog().targets.length*realDatasets().length<=10000&&realDatasets().length<=32;
  const fmt=value=>formatStatistic(value,formatAbundance,'raw',$('number-format').value);
  const requestKey=(prefix,extra=[])=>context.controlKey(['threshold','percentiles','quantiles',...extra].map(suffix=>prefix+'-'+suffix))+'|'+context.normalization()+'|'+context.metric();
  const overviewKey=()=>requestKey('overview');
  const datasetsKey=()=>requestKey('dataset-list',['scope','targets']);
  const focusKey=()=>requestKey('target-focus',['scope','limit','name']);
  function updateActions() {
    if (!catalog()) return;
    context.actionState('compute-overview',overview?.key!==overviewKey());
    let ids;
    try {ids=orderedDatasets();} catch {ids=datasetScope();}
    const same=datasetResults?.key===datasetsKey();
    context.actionState('compute-dataset-page',!same||ids.slice(datasetPage*PAGE,(datasetPage+1)*PAGE).some(id=>!datasetResults.rows.has(id)));
    context.actionState('compute-datasets-all',!same||ids.some(id=>!datasetResults.rows.has(id)));
    context.actionState('load-target-focus',focus?.key!==focusKey());
  }
  function options(prefix) {
    const text=$(prefix+'-threshold').value.trim(),isInteger=integer();
    if (!(isInteger?/^\d+$/:/^\d+(\.\d+)?$/).test(text)) throw new Error('Threshold must be a nonnegative '+(isInteger?'integer.':'number.'));
    const threshold=isInteger?BigInt(text):Number(text);
    if (!isInteger&&!Number.isFinite(threshold)) throw new Error('Threshold must be finite.');
    return {integer:isInteger,threshold,percentiles:parsePercentiles($(prefix+'-percentiles').value),method:$(prefix+'-quantiles').value,
      sampleSize:4096,normalization:context.normalization(),metric:context.metric()};
  }
  const distribution=opts=>new Distribution(opts);
  function validate(targets,datasets,opts) {
    if (!targets.length||!datasets.length) throw new Error('Choose at least one target and one real dataset.');
    for (const id of datasets) normalize(opts.integer?0n:0,catalog().targets[targets[0]],catalog().datasets[id],opts.normalization);
    if (['coverage','length-cpm'].includes(opts.normalization)) for (const id of targets)
      normalize(0n,catalog().targets[id],catalog().datasets[datasets[0]],opts.normalization);
  }
  async function scanDataset(dataset,targets,opts,version,visit) {
    for (let start=0;start<targets.length;start+=4096) {
      context.checkVersion(version);
      const ids=targets.slice(start,start+4096),values=await context.read(dataset,ids);
      context.checkVersion(version);
      if (values.length!==ids.length) throw new Error('Incomplete abundance values.');
      values.forEach((text,i)=>{
        const raw=context.parseValue(text),value=raw===null?null:normalize(raw,catalog().targets[ids[i]],catalog().datasets[dataset],opts.normalization);
        visit(value,raw,ids[i]);
      });
    }
  }
  async function run(action) {
    if (!catalog()||context.busy()) return;
    context.setBusy(true);
    try {await action(context.version());return true;}
    catch(error) {context.status(error.message,error.message!=='Operation canceled.');return false;}
    finally {context.setBusy(false);render();}
  }
  function cards(root,fields) {
    $(root).replaceChildren(...fields.map(([name,value,key])=>{
      const card=document.createElement('div'),label=document.createElement('small'),number=document.createElement('strong');
      if (key) card.dataset.statistic=key;
      label.textContent=name;number.textContent=typeof value==='string'&&value.endsWith('%')?value:fmt(value);number.title=value?.toString()??'Not computed / missing';
      card.append(label,number);return card;
    }));
  }
  function statisticFields(stats,opts,unit) {
    return [['Total expression',stats?.total,'total'],['Known '+unit,stats?.known,'known'],['Missing '+unit,stats?.missing,'missing'],
      ['Nonzero '+unit,stats?.nonzero,'nonzero'],['Mean',stats?.mean,'mean'],['Median',stats?.median,'median'],
      ['Above '+(opts?.threshold??0),stats?.above,'above'],['% above threshold',stats?.above_percent==null?null:percentage(stats.above_percent),'above_percent'],
      ...((opts?.percentiles)||[10,90,95,99]).map(p=>['P'+p,stats?.percentiles[p],'p'+p])];
  }
  function renderOverview() {
    if (!catalog()) return;
    const datasets=realDatasets(),cells=BigInt(catalog().targets.length)*BigInt(datasets.length);
    const fields=[[context.targetKind()+'s',catalog().targets.length,'targets'],['Real datasets',datasets.length,'datasets'],['Matrix cells',cells,'cells'],['Index bytes',context.file().size,'bytes']];
    cards('overview-cards',fields.concat(statisticFields(overview?.statistics,overview?.options,'cells'),[['Zero cells',overview?overview.statistics.known-overview.statistics.nonzero:null,'zero']]));
    if (overview) {
      const stats=overview.statistics;
      $('overview-details').textContent=[stats.quantiles_approximate?'Approximate quantiles':'',catalog().coverage?.complete===false?'Incomplete dataset coverage':''].filter(Boolean).join(' · ');
    } else $('overview-details').textContent='';
    $('download-overview').disabled=context.busy()||!overview;
    updateActions();
  }
  async function computeOverview() {
    if (overview?.key===overviewKey()) return true;
    return await run(async version=>{
      const opts=options('overview'),targets=targetIds(),datasets=realDatasets();validate(targets,datasets,opts);
      const dist=distribution(opts);
      context.progress(datasets.length*Math.ceil(targets.length/4096),'Global summary');
      for (const [i,id] of datasets.entries()) {
        context.status(`Global overview: dataset ${i+1} of ${datasets.length.toLocaleString()} · ${catalog().datasets[id].name}`);
        await scanDataset(id,targets,opts,version,value=>dist.add(value));
      }
      context.checkVersion(version);overview={statistics:dist.result(opts.percentiles),options:opts,key:overviewKey()};
      context.status('');
    });
  }
  function datasetScope() {return realDatasets($('dataset-list-scope').value);}
  function orderedDatasets() {
    const ids=datasetScope(),field=$('dataset-list-sort').value,direction=$('dataset-list-direction').value==='asc'?1:-1;
    if (!['reference','name'].includes(field)&&ids.some(id=>!datasetResults?.rows.has(id))) throw new Error('Compute all datasets in scope before sorting by expression statistics.');
    return ids.sort((a,b)=>{
      if (field==='reference') return a-b;
      if (field==='name') return direction*catalog().datasets[a].name.localeCompare(catalog().datasets[b].name,undefined,{numeric:true})||a-b;
      const x=datasetResults.rows.get(a)[field],y=datasetResults.rows.get(b)[field];
      if (x===null||y===null) return x===null?(y===null?a-b:1):-1;
      return direction*compareStatistics(x,y)||a-b;
    });
  }
  function datasetColumns() {
    const opts=datasetResults?.options;
    return [['known','Known targets'],['missing','Missing'],['total','Total'],['nonzero','Expressed targets'],['above','Targets > X'],['mean','Mean'],['median','Median'],['min','Min'],['max','Max'],['stddev','Std dev'],...(opts?.percentiles||[10,90,95,99]).map(p=>['p'+p,'P'+p])];
  }
  const fieldValue=(stats,field)=>field.startsWith('p')?stats?.percentiles[Number(field.slice(1))]:stats?.[field];
  function table(root,labels,rows) {
    const table=document.createElement('table'),head=document.createElement('thead'),header=document.createElement('tr');
    for (const label of labels) {const th=document.createElement('th');th.textContent=label;header.append(th);}
    head.append(header);table.append(head);const body=document.createElement('tbody');
    for (const row of rows) body.append(row);table.append(body);$(root).replaceChildren(table);
  }
  function linkCell(name,action) {
    const th=document.createElement('th'),button=document.createElement('button');button.className='table-link';button.textContent=name;
    button.disabled=context.busy();button.addEventListener('click',action);th.append(button);return th;
  }
  function renderDatasets() {
    if (!catalog()) return;
    let ids;
    try {ids=orderedDatasets();} catch(error) {context.status(error.message,true);return;}
    datasetPage=Math.min(datasetPage,Math.max(0,Math.ceil(ids.length/PAGE)-1));
    const page=ids.slice(datasetPage*PAGE,(datasetPage+1)*PAGE),columns=datasetColumns();
    table('datasets-table',['Dataset',...columns.map(([_,name])=>name)],page.map(id=>{
      const row=document.createElement('tr');row.dataset.dataset=id;
      row.append(linkCell(catalog().datasets[id].name,()=>context.focusDataset(id)));
      for (const [field] of columns) {const td=document.createElement('td');td.dataset.statistic=field;td.textContent=fmt(fieldValue(datasetResults?.rows.get(id),field));row.append(td);}
      return row;
    }));
    const computed=ids.filter(id=>datasetResults?.rows.has(id)).length;
    $('datasets-details').textContent=`${computed.toLocaleString()} / ${ids.length.toLocaleString()} datasets calculated`;
    $('datasets-page-label').textContent=`Datasets ${ids.length?datasetPage*PAGE+1:0}–${Math.min((datasetPage+1)*PAGE,ids.length)} of ${ids.length.toLocaleString()}`;
    $('datasets-prev').disabled=context.busy()||datasetPage===0;$('datasets-next').disabled=context.busy()||(datasetPage+1)*PAGE>=ids.length;
    $('download-datasets').disabled=context.busy()||!computed;
    updateActions();
  }
  async function computeDatasets(all=false) {
    if (datasetResults?.key===datasetsKey()) {
      const ids=orderedDatasets(),requested=all?ids:ids.slice(datasetPage*PAGE,(datasetPage+1)*PAGE);
      if (requested.every(id=>datasetResults.rows.has(id))) {renderDatasets();return;}
    }
    await run(async version=>{
      const opts=options('dataset-list'),targets=targetIds($('dataset-list-targets').value),ids=datasetScope();validate(targets,ids,opts);
      if (datasetResults?.key!==datasetsKey()) datasetResults={rows:new Map(),options:opts,key:datasetsKey()};
      // Statistical sort is enabled only after full scope coverage; scan reference order otherwise.
      const ordered=['reference','name'].includes($('dataset-list-sort').value)?orderedDatasets():ids;
      const requested=all?ordered:ordered.slice(datasetPage*PAGE,(datasetPage+1)*PAGE);
      context.progress(requested.filter(id=>!datasetResults.rows.has(id)).length*Math.ceil(targets.length/4096),'Dataset list');
      for (const [i,id] of requested.entries()) {
        if (datasetResults.rows.has(id)) continue;
        context.status(`Dataset list: ${i+1} of ${requested.length} · ${catalog().datasets[id].name}`);
        const dist=distribution(opts);await scanDataset(id,targets,opts,version,value=>dist.add(value));context.checkVersion(version);
        datasetResults.rows.set(id,dist.result(opts.percentiles));
      }
      context.status('');
    });
  }
  function targetId(text) {
    if (/^#[1-9]\d*$/.test(text)) {const id=Number(text.slice(1))-1;if (Number.isSafeInteger(id)&&id<catalog().targets.length) return id;}
    const id=catalog().targets.findIndex(target=>target.name===text);
    if (id<0) throw new Error('Unknown target: '+text);return id;
  }
  function renderFocus() {
    if (!catalog()) return;
    $('download-target-focus').disabled=context.busy()||!focus;$('download-target-svg').disabled=context.busy()||!focus;
    $('download-target-page').disabled=context.busy()||!focus?.page;
    $('target-focus-prev').disabled=context.busy()||!focus||focusPage===0;
    $('target-focus-next').disabled=context.busy()||!focus||(focusPage+1)*PAGE>=focus.datasets.length;
    updateActions();
    if (!focus) {$('target-focus-title').textContent=context.targetKind()+' focus';cards('target-focus-cards',[]);$('target-focus-chart').hidden=true;$('target-focus-table').replaceChildren();$('target-focus-page-label').textContent='';return;}
    const target=catalog().targets[focus.id],stats=focus.statistics;
    $('target-focus-title').textContent=context.targetKind()+' focus · '+target.name+' (#'+(focus.id+1)+')';
    cards('target-focus-cards',statisticFields(stats,focus.options,'datasets'));
    const rows=focus.ranking.rows.map(row=>({...row,name:catalog().datasets[row.id].name,percent:share(row.value,focus.ranking.total)}));
    const svg=barplot({rows,title:'Most expressive datasets · '+target.name,subtitle:'% of target expression',format:fmt,valueLabel:context.metricLabel(),denominatorLabel:'target signal',accent:context.colors().main});
    $('target-focus-chart').replaceChildren(svg);$('target-focus-chart').hidden=false;
    table('target-focus-table',['Dataset','Expression','Fold change','Log₂ fold change'],(focus.page||[]).map(entry=>{
      const row=document.createElement('tr');row.append(linkCell(catalog().datasets[entry.dataset].name,()=>context.focusDataset(entry.dataset)));
      for (const [field,value] of [['value',entry.value],['fold',foldChange(entry.value,stats.total,stats.known)],['log-fold',foldChange(entry.value,stats.total,stats.known,'log-fold')]]) {
        const td=document.createElement('td');td.dataset.statistic=field;td.textContent=fmt(value);td.title=entry.raw===null?'Missing estimate':'Raw value: '+entry.raw;row.append(td);
      }return row;
    }));
    $('target-focus-details').textContent=stats.quantiles_approximate?'Approximate quantiles':'';
    $('target-focus-page-label').textContent=`Datasets ${focusPage*PAGE+1}–${Math.min((focusPage+1)*PAGE,focus.datasets.length)} of ${focus.datasets.length.toLocaleString()}`;
  }
  async function loadFocus() {
    if (focus?.key===focusKey()) return;
    await run(async version=>{
      const id=targetId($('target-focus-name').value.trim()),opts=options('target-focus'),datasets=realDatasets($('target-focus-scope').value);validate([id],datasets,opts);
      const ranking=new Ranking(Number($('target-focus-limit').value),opts.integer),dist=distribution(opts),page=[];
      context.progress(datasets.length,'Target profile');
      for (const [i,dataset] of datasets.entries()) {
        context.status(`Target focus: dataset ${i+1} of ${datasets.length.toLocaleString()} · ${catalog().datasets[dataset].name}`);
        await scanDataset(dataset,[id],opts,version,(value,raw)=>{dist.add(value);ranking.add(dataset,value);if (i<PAGE) page.push({dataset,value,raw});});
      }
      context.checkVersion(version);focusPage=0;focus={id,datasets,options:opts,statistics:dist.result(opts.percentiles),ranking:ranking.result(),page,key:focusKey()};
      context.status('');
    });
  }
  async function focusPageLoad(direction) {
    if (!focus||context.busy()) return;
    await run(async version=>{
      const next=focusPage+direction,page=[];
      if (next<0||next*PAGE>=focus.datasets.length) return;
      context.progress(Math.min(PAGE,focus.datasets.length-next*PAGE),'Target profile');
      for (const dataset of focus.datasets.slice(next*PAGE,(next+1)*PAGE)) await scanDataset(dataset,[focus.id],focus.options,version,(value,raw)=>page.push({dataset,value,raw}));
      context.checkVersion(version);focusPage=next;focus.page=page;
    });
  }
  const csvCell=value=>'"'+String(value??'').replaceAll('"','""')+'"';
  function download(rows,name,type='text/csv;charset=utf-8') {
    const text=typeof rows==='string'?rows:rows.map(row=>row.map(csvCell).join(',')).join('\r\n')+'\r\n';
    const url=URL.createObjectURL(new Blob([text],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),10000);
  }
  function statsCsv(name,stats,opts,identity) {
    const fields=['total','known','missing','nonzero','above','above_percent','mean','median','min','max','stddev'];
    download([['name','metric','normalization','threshold','quantiles_approximate','sample_n',...fields,...opts.percentiles.map(p=>'p'+p)],
      [identity,opts.metric,opts.normalization,opts.threshold,stats.quantiles_approximate,stats.sample_n,...fields.map(field=>stats[field]),...opts.percentiles.map(p=>stats.percentiles[p])]],name);
  }
  function render() {renderOverview();renderDatasets();renderFocus();}
  function selectionChanged() {datasetResults=undefined;datasetPage=0;$('dataset-list-sort').value='reference';focus=undefined;focusPage=0;render();}
  function resetSource() {overview=undefined;$('dataset-list-sort').value='reference';selectionChanged();}
  $('compute-overview').addEventListener('click',computeOverview);
  $('compute-dataset-page').addEventListener('click',()=>computeDatasets());$('compute-datasets-all').addEventListener('click',()=>computeDatasets(true));
  $('load-target-focus').addEventListener('click',loadFocus);
  for (const [id,direction] of [['datasets-prev',-1],['datasets-next',1]]) $(id).addEventListener('click',()=>{if (!context.busy()) {datasetPage+=direction;renderDatasets();computeDatasets();}});
  for (const [id,direction] of [['target-focus-prev',-1],['target-focus-next',1]]) $(id).addEventListener('click',()=>focusPageLoad(direction));
  for (const suffix of ['threshold','percentiles','quantiles']) $('overview-'+suffix).addEventListener('change',()=>{if (overview?.key!==overviewKey()) overview=undefined;renderOverview();});
  for (const suffix of ['scope','targets','threshold','percentiles','quantiles']) $('dataset-list-'+suffix).addEventListener('change',()=>{if (datasetResults?.key===datasetsKey()) return;datasetResults=undefined;datasetPage=0;$('dataset-list-sort').value='reference';renderDatasets();});
  for (const suffix of ['sort','direction']) $('dataset-list-'+suffix).addEventListener('change',()=>{
    try {orderedDatasets();} catch(error) {context.status(error.message,true);$('dataset-list-sort').value='reference';}
    datasetPage=0;renderDatasets();
  });
  for (const suffix of ['scope','threshold','percentiles','quantiles','limit','name']) $('target-focus-'+suffix).addEventListener('change',()=>{if (focus?.key!==focusKey()) focus=undefined;renderFocus();});
  $('download-overview').addEventListener('click',()=>{if (overview&&!context.busy()) statsCsv('expresso-matrix-summary.csv',overview.statistics,overview.options,'whole matrix');});
  $('download-datasets').addEventListener('click',()=>{
    if (!datasetResults||context.busy()) return;
    const cols=datasetColumns(),ids=orderedDatasets(),computed=ids.filter(id=>datasetResults.rows.has(id)),opts=datasetResults.options;
    download([['dataset','metric','normalization','target_scope','dataset_scope','scope_datasets','calculated_datasets','threshold','quantiles_approximate',...cols.map(([field])=>field)],
      ...computed.map(id=>{const stats=datasetResults.rows.get(id);return [catalog().datasets[id].name,opts.metric,opts.normalization,$('dataset-list-targets').value,$('dataset-list-scope').value,ids.length,computed.length,opts.threshold,stats.quantiles_approximate,...cols.map(([field])=>fieldValue(stats,field))];})],'expresso-dataset-list.csv');
  });
  $('download-target-focus').addEventListener('click',()=>{if (focus&&!context.busy()) statsCsv('expresso-target-summary.csv',focus.statistics,focus.options,catalog().targets[focus.id].name);});
  $('download-target-svg').addEventListener('click',()=>{if (focus&&!context.busy()) download(new XMLSerializer().serializeToString($('target-focus-chart').querySelector('svg')),'expresso-target-profile.svg','image/svg+xml;charset=utf-8');});
  $('download-target-page').addEventListener('click',()=>{
    if (!focus||context.busy()) return;
    download([['target_id','target_name','dataset','metric','normalization','raw_value','expression','fold_change','log2_fold_change'],
      ...focus.page.map(entry=>[focus.id+1,catalog().targets[focus.id].name,catalog().datasets[entry.dataset].name,focus.options.metric,focus.options.normalization,entry.raw,entry.value,
        foldChange(entry.value,focus.statistics.total,focus.statistics.known),foldChange(entry.value,focus.statistics.total,focus.statistics.known,'log-fold')])],'expresso-target-visible-profile.csv');
  });
  return {small,render,renderOverview,resetSource,selectionChanged,computeOverview,computeDatasets,loadFocus,updateActions,
    focusTarget:async id=>{$('target-focus-name').value='#'+(id+1);context.activate('target-focus');await loadFocus();},
    enter:async name=>{
      if (name==='overview') renderOverview();
      if (name==='datasets') {renderDatasets();await computeDatasets();}
      if (name==='target-focus') {renderFocus();if (!focus&&small()) await loadFocus();}
    }};
}
