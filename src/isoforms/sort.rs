//! Bounded-memory, bounded-fan-in external sorting of reference occurrences.
use anyhow::{Context, Result};
use rayon::prelude::*;
use std::{
    cmp::Reverse,
    collections::BinaryHeap,
    fs::{self, File},
    io::{BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct Event {
    pub kmer: u64,
    // Even: exon ID. Odd: transcript ID. IDs are zero-based.
    pub member: u64,
    pub count: u64,
}
fn write(out: &mut impl Write, e: Event) -> Result<()> {
    for value in [e.kmer, e.member, e.count] {
        out.write_all(&value.to_le_bytes())?;
    }
    Ok(())
}
fn read(input: &mut impl Read) -> Result<Option<Event>> {
    let mut bytes = [0; 24];
    if input.read(&mut bytes[..1])? == 0 {
        return Ok(None);
    }
    input.read_exact(&mut bytes[1..])?;
    Ok(Some(Event {
        kmer: u64::from_le_bytes(bytes[..8].try_into()?),
        member: u64::from_le_bytes(bytes[8..16].try_into()?),
        count: u64::from_le_bytes(bytes[16..].try_into()?),
    }))
}

pub(super) struct Sorter<'a> {
    root: PathBuf,
    capacity: usize,
    pending: Vec<Event>,
    files: Vec<PathBuf>,
    serial: usize,
    pool: &'a rayon::ThreadPool,
}
impl<'a> Sorter<'a> {
    pub fn new(root: &Path, capacity: usize, pool: &'a rayon::ThreadPool) -> Result<Self> {
        fs::create_dir_all(root)?;
        Ok(Self {
            root: root.into(),
            capacity: capacity.max(1),
            pending: Vec::new(),
            files: Vec::new(),
            serial: 0,
            pool,
        })
    }
    fn path(&mut self) -> PathBuf {
        self.serial += 1;
        self.root.join(format!("{}.bin", self.serial))
    }
    pub fn root(&self) -> &Path {
        &self.root
    }
    pub fn capacity(&self) -> usize {
        self.capacity
    }
    pub fn add_runs(&mut self, files: Vec<PathBuf>) {
        self.files.extend(files);
    }
    pub fn into_runs(mut self) -> Result<Vec<PathBuf>> {
        self.flush()?;
        Ok(self.files)
    }
    pub fn push(&mut self, kmer: u64, member: u64) -> Result<()> {
        self.pending.push(Event {
            kmer,
            member,
            count: 1,
        });
        if self.pending.len() >= self.capacity {
            self.flush()?;
        }
        Ok(())
    }
    fn flush(&mut self) -> Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        self.pool.install(|| self.pending.par_sort_unstable());
        let path = self.path();
        let mut out = BufWriter::new(File::create(&path)?);
        let mut last: Option<Event> = None;
        for &e in &self.pending {
            if let Some(old) = last.as_mut()
                && (old.kmer, old.member) == (e.kmer, e.member)
            {
                old.count = old
                    .count
                    .checked_add(e.count)
                    .context("reference multiplicity overflow")?;
            } else {
                if let Some(old) = last {
                    write(&mut out, old)?;
                }
                last = Some(e);
            }
        }
        if let Some(old) = last {
            write(&mut out, old)?;
        }
        out.flush()?;
        self.pending.clear();
        self.files.push(path);
        Ok(())
    }
    pub fn finish(self, visit: impl FnMut(Event) -> Result<()>) -> Result<()> {
        self.finish_with_fan_in(64, visit)
    }
    fn finish_with_fan_in(
        mut self,
        fan_in: usize,
        visit: impl FnMut(Event) -> Result<()>,
    ) -> Result<()> {
        self.flush()?;
        drop(std::mem::take(&mut self.pending));
        while self.files.len() > fan_in {
            let files = std::mem::take(&mut self.files);
            for group in files.chunks(fan_in) {
                let path = self.path();
                let mut out = BufWriter::new(File::create(&path)?);
                merge(group, |e| write(&mut out, e))?;
                out.flush()?;
                for input in group {
                    fs::remove_file(input)?;
                }
                self.files.push(path);
            }
        }
        merge(&self.files, visit)
    }
}
fn merge(files: &[PathBuf], mut visit: impl FnMut(Event) -> Result<()>) -> Result<()> {
    let mut readers: Vec<_> = files
        .iter()
        .map(|path| Ok(BufReader::with_capacity(65536, File::open(path)?)))
        .collect::<Result<_>>()?;
    let mut heap = BinaryHeap::new();
    for (id, reader) in readers.iter_mut().enumerate() {
        if let Some(e) = read(reader)? {
            heap.push(Reverse((e, id)));
        }
    }
    let mut last: Option<Event> = None;
    while let Some(Reverse((e, id))) = heap.pop() {
        if let Some(old) = last.as_mut()
            && (old.kmer, old.member) == (e.kmer, e.member)
        {
            old.count = old
                .count
                .checked_add(e.count)
                .context("reference multiplicity overflow")?;
        } else {
            if let Some(old) = last {
                visit(old)?;
            }
            last = Some(e);
        }
        if let Some(e) = read(&mut readers[id])? {
            heap.push(Reverse((e, id)));
        }
    }
    if let Some(old) = last {
        visit(old)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn spills_and_multilevel_merges_aggregate_multiplicity() {
        let dir = tempfile::tempdir().unwrap();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        let mut sorter = Sorter::new(dir.path(), 2, &pool).unwrap();
        for kmer in (0..25).rev() {
            for _ in 0..3 {
                sorter.push(kmer, kmer % 5).unwrap();
            }
        }
        let mut found = Vec::new();
        sorter
            .finish_with_fan_in(3, |e| {
                found.push(e);
                Ok(())
            })
            .unwrap();
        assert_eq!(
            found,
            (0..25)
                .map(|kmer| Event {
                    kmer,
                    member: kmer % 5,
                    count: 3
                })
                .collect::<Vec<_>>()
        );
    }
}
