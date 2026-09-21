use anyhow::Result;
use clap::ValueEnum;
use serde::Serialize;
use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    path::Path,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Format {
    Compact,
    Csv,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Compression {
    None,
    Gz,
    Xz,
    Zstd,
}
impl Compression {
    pub fn vector_suffix(self) -> &'static str {
        match self {
            Self::None => ".eab",
            Self::Gz => ".eab.gz",
            Self::Xz => ".eab.xz",
            Self::Zstd => ".eab.zst",
        }
    }
    pub fn suffix(self) -> &'static str {
        match self {
            Self::None => ".csv",
            Self::Gz => ".csv.gz",
            Self::Xz => ".csv.xz",
            Self::Zstd => ".csv.zst",
        }
    }
}

enum Encoder {
    Plain(BufWriter<File>),
    Gz(flate2::write::GzEncoder<BufWriter<File>>),
    Xz(xz2::write::XzEncoder<BufWriter<File>>),
    Zstd(zstd::stream::write::Encoder<'static, BufWriter<File>>),
}
impl Write for Encoder {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(buf),
            Self::Gz(w) => w.write(buf),
            Self::Xz(w) => w.write(buf),
            Self::Zstd(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Gz(w) => w.flush(),
            Self::Xz(w) => w.flush(),
            Self::Zstd(w) => w.flush(),
        }
    }
}

pub fn csv(
    path: &Path,
    compression: Compression,
    write: impl FnOnce(&mut csv::Writer<&mut dyn Write>) -> Result<()>,
) -> Result<()> {
    encoded(path, compression, |encoder| {
        let mut csv = csv::Writer::from_writer(encoder);
        write(&mut csv)?;
        csv.flush()?;
        Ok(())
    })
}

pub fn encoded(
    path: &Path,
    compression: Compression,
    write: impl FnOnce(&mut dyn Write) -> Result<()>,
) -> Result<()> {
    let file = BufWriter::with_capacity(1 << 20, File::create(path)?);
    let mut encoder = match compression {
        Compression::None => Encoder::Plain(file),
        Compression::Gz => Encoder::Gz(flate2::write::GzEncoder::new(
            file,
            flate2::Compression::fast(),
        )),
        Compression::Xz => Encoder::Xz(xz2::write::XzEncoder::new(file, 3)),
        Compression::Zstd => Encoder::Zstd(zstd::stream::write::Encoder::new(file, 3)?),
    };
    write(&mut encoder)?;
    let mut file = match encoder {
        Encoder::Plain(file) => file,
        Encoder::Gz(w) => w.finish()?,
        Encoder::Xz(w) => w.finish()?,
        Encoder::Zstd(w) => w.finish()?,
    };
    file.flush()?;
    Ok(())
}

pub fn reader(path: &Path) -> Result<Box<dyn Read>> {
    let mut input = BufReader::new(File::open(path)?);
    let magic = input.fill_buf()?;
    Ok(if magic.starts_with(b"\x28\xb5\x2f\xfd") {
        Box::new(zstd::stream::read::Decoder::new(input)?)
    } else if magic.starts_with(b"\x1f\x8b") {
        Box::new(flate2::read::MultiGzDecoder::new(input))
    } else if magic.starts_with(b"\xfd7zXZ\0") {
        Box::new(xz2::read::XzDecoder::new_multi_decoder(input))
    } else {
        Box::new(input)
    })
}
