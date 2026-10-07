//! Minimal localhost server for the bundled, entirely offline WebAssembly app.
use anyhow::Result;
use std::{
    io::{Read, Write},
    net::{TcpListener, TcpStream},
};

fn request(mut stream: TcpStream) -> Result<()> {
    stream.set_read_timeout(Some(std::time::Duration::from_secs(5)))?;
    stream.set_write_timeout(Some(std::time::Duration::from_secs(5)))?;
    let mut bytes = Vec::new();
    let mut chunk = [0; 1024];
    while !bytes.windows(4).any(|s| s == b"\r\n\r\n") && bytes.len() < 8192 {
        let n = stream.read(&mut chunk)?;
        if n == 0 {
            return Ok(());
        }
        bytes.extend_from_slice(&chunk[..n]);
    }
    let request = std::str::from_utf8(&bytes)?.lines().next().unwrap_or("");
    let mut parts = request.split_ascii_whitespace();
    let method = parts.next().unwrap_or("");
    let path = parts.next().unwrap_or("").split('?').next().unwrap_or("");
    let asset: Option<(&str, &[u8])> = match path {
        "/" | "/index.html" => Some((
            "text/html; charset=utf-8",
            include_bytes!("../viewer/index.html"),
        )),
        "/style.css" => Some((
            "text/css; charset=utf-8",
            include_bytes!("../viewer/style.css"),
        )),
        "/app.js" => Some((
            "text/javascript; charset=utf-8",
            include_bytes!("../viewer/app.js"),
        )),
        "/worker.js" => Some((
            "text/javascript; charset=utf-8",
            include_bytes!("../viewer/worker.js"),
        )),
        "/display.js" => Some((
            "text/javascript; charset=utf-8",
            include_bytes!("../viewer/display.js"),
        )),
        "/expresso_viewer.wasm" => Some((
            "application/wasm",
            include_bytes!("../viewer/expresso_viewer.wasm"),
        )),
        _ => None,
    };
    let (status, mime, body) = if !["GET", "HEAD"].contains(&method) {
        (
            "405 Method Not Allowed",
            "text/plain",
            b"Method not allowed".as_slice(),
        )
    } else if let Some((mime, body)) = asset {
        ("200 OK", mime, body)
    } else {
        ("404 Not Found", "text/plain", b"Not found".as_slice())
    };
    write!(
        stream,
        "HTTP/1.1 {status}\r\nContent-Type: {mime}\r\nContent-Length: {}\r\nConnection: close\r\nCache-Control: no-store\r\n\r\n",
        body.len()
    )?;
    if method != "HEAD" {
        stream.write_all(body)?;
    }
    Ok(())
}

pub fn serve(port: u16) -> Result<()> {
    let listener = TcpListener::bind((std::net::Ipv4Addr::LOCALHOST, port))?;
    eprintln!(
        "EXPRESSO viewer: http://{} (Ctrl+C to stop)",
        listener.local_addr()?
    );
    for stream in listener.incoming() {
        let stream = stream?;
        std::thread::spawn(move || {
            if let Err(error) = request(stream) {
                eprintln!("Viewer request: {error}");
            }
        });
    }
    Ok(())
}
