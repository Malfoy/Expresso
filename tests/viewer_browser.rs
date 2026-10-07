//! Optional real-browser regression: run a WebDriver server on localhost:4444 first.
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::{BufRead, BufReader, Read, Write},
    net::TcpStream,
    process::{Child, Command, Stdio},
    thread,
    time::Duration,
};

fn webdriver(method: &str, path: &str, value: Value) -> Value {
    let bytes = serde_json::to_vec(&value).unwrap();
    let address = std::env::var("EXPRESSO_WEBDRIVER_ADDRESS")
        .unwrap_or_else(|_| "127.0.0.1:4444".to_string());
    let mut stream =
        TcpStream::connect(&address).expect("Start geckodriver --host 127.0.0.1 --port 4444 first");
    stream
        .set_read_timeout(Some(Duration::from_secs(60)))
        .unwrap();
    write!(stream,"{method} {path} HTTP/1.1\r\nHost: {address}\r\nConnection: close\r\nContent-Type: application/json\r\nContent-Length: {}\r\n\r\n",bytes.len()).unwrap();
    stream.write_all(&bytes).unwrap();
    let mut response = String::new();
    stream.read_to_string(&mut response).unwrap();
    let (_, body) = response.split_once("\r\n\r\n").expect("HTTP response");
    let value: Value = serde_json::from_str(body).expect(body);
    assert!(value["value"]["error"].is_null(), "{value}");
    value["value"].clone()
}
struct Browser {
    path: String,
    server: Child,
}
impl Browser {
    fn js(&self, script: &str) -> Value {
        webdriver(
            "POST",
            &format!("{}/execute/sync", self.path),
            json!({"script":script,"args":[]}),
        )
    }
    fn click(&self, selector: &str) {
        self.js(&format!(
            "document.querySelector({}).click()",
            serde_json::to_string(selector).unwrap()
        ));
    }
    fn wait(&self, condition: &str) {
        for _ in 0..300 {
            if self.js(condition) == true {
                return;
            }
            thread::sleep(Duration::from_millis(100));
        }
        panic!(
            "Timed out: {condition}, status {}",
            self.js("return document.getElementById('status').textContent")
        );
    }
}
impl Drop for Browser {
    fn drop(&mut self) {
        let _ = self.server.kill();
        let _ = self.server.wait();
        // Avoid assertions during unwinding if the external driver stopped.
        let _ = std::panic::catch_unwind(|| webdriver("DELETE", &self.path, Value::Null));
    }
}

#[test]
#[ignore = "requires Firefox and geckodriver listening on localhost:4444"]
fn large_selection_reads_only_visible_counts_and_streams_full_csv() {
    // Use the project directory so sandboxed Firefox can access the selected local file.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir_in(root.join("target")).unwrap();
    let dir = temp.path();
    let names: Vec<_> = (1..=300).map(|i| format!("g{i:03}")).collect();
    let mut hash = Sha256::new();
    hash.update(b"EXPRESSO-reference-v1\0");
    hash.update(300u64.to_le_bytes());
    let mut table = String::from("gene_id,gene_name,length,unique_kmers\n");
    for (i, name) in names.iter().enumerate() {
        hash.update((name.len() as u64).to_le_bytes());
        hash.update(name);
        hash.update(7u64.to_le_bytes());
        hash.update(1u64.to_le_bytes());
        table.push_str(&format!("{},{name},7,1\n", i + 1));
    }
    let reference: [u8; 32] = hash.finalize().into();
    let hex: String = reference.iter().map(|b| format!("{b:02x}")).collect();
    let mut raw = b"EXPRAB01".to_vec();
    raw.extend([2, 0, 0, 0]);
    raw.extend(300u64.to_le_bytes());
    raw.extend(4u32.to_le_bytes());
    raw.extend(75u64.to_le_bytes());
    raw.extend(u64::MAX.to_le_bytes());
    raw.extend(reference);
    for value in [0u64, 1, 2, u64::MAX] {
        raw.extend(value.to_le_bytes());
    }
    raw.extend([0xe4; 75]);
    raw.extend(crc32fast::hash(&raw).to_le_bytes());
    fs::write(dir.join("reference.csv"), table).unwrap();
    fs::write(dir.join("vector.eab"), raw).unwrap();
    // g001 has signal only beyond the first 25-dataset page. Global zero-row
    // filtering must keep it; sample_026 also breaks otherwise identical totals.
    let mut alternate = fs::read(dir.join("vector.eab")).unwrap();
    alternate[104] = 0xe5;
    let payload_end = alternate.len() - 4;
    let checksum = crc32fast::hash(&alternate[..payload_end]);
    alternate[payload_end..].copy_from_slice(&checksum.to_le_bytes());
    fs::write(dir.join("alternate.eab"), alternate).unwrap();
    let manifest = json!({"format":"compact","compact_format_version":1,"level":"gene","reference":{"file":"reference.csv","sha256":hex,"targets":300},
        "datasets":(1..=150).map(|i|json!({"name":format!("sample_{i:03}"),"output":if i==26 {"alternate.eab"} else {"vector.eab"}})).collect::<Vec<_>>(),"global":{"output":"vector.eab"},"coverage":{"complete":true}});
    fs::write(
        dir.join("manifest.json"),
        serde_json::to_vec(&manifest).unwrap(),
    )
    .unwrap();
    let index = dir.join("large.eai");
    assert!(
        Command::new(env!("CARGO_BIN_EXE_expresso"))
            .args(["pack", "--input"])
            .arg(dir)
            .arg("--output")
            .arg(&index)
            .output()
            .unwrap()
            .status
            .success()
    );
    let mut server = Command::new(env!("CARGO_BIN_EXE_expresso"))
        .args(["view", "--port", "0"])
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut line = String::new();
    BufReader::new(server.stderr.take().unwrap())
        .read_line(&mut line)
        .unwrap();
    let url = line.split_ascii_whitespace().nth(2).unwrap();
    let session = webdriver(
        "POST",
        "/session",
        json!({"capabilities":{"alwaysMatch":{"browserName":"firefox","moz:firefoxOptions":{"args":["-headless"]}}}}),
    );
    let browser = Browser {
        path: format!("/session/{}", session["sessionId"].as_str().unwrap()),
        server,
    };
    webdriver("POST", &format!("{}/url", browser.path), json!({"url":url}));
    browser.js("window.queries=[];const send=Worker.prototype.postMessage;Worker.prototype.postMessage=function(message){window.queries.push({offset:message.offset,targets:message.request.selected.length});return send.call(this,message);}");
    let input = webdriver(
        "POST",
        &format!("{}/element", browser.path),
        json!({"using":"css selector","value":"#file"}),
    );
    let element = input["element-6066-11e4-a52e-4f735466cecf"]
        .as_str()
        .unwrap();
    webdriver(
        "POST",
        &format!("{}/element/{element}/value", browser.path),
        json!({"text":index.canonicalize().unwrap().to_str().unwrap()}),
    );
    browser.wait("return !document.getElementById('controls').hidden");
    assert_eq!(
        browser.js("return window.queries.length"),
        0,
        "opening the file must not read any abundance vectors"
    );
    browser.click("#select-targets");
    browser.click("#select-datasets");
    browser.click("#draw");
    browser.wait("return !document.getElementById('download').disabled");
    assert_eq!(
        browser.js("return document.querySelectorAll('#heatmap td').length"),
        2500
    );
    assert_eq!(browser.js("return window.queries.length"), 25);
    assert_eq!(
        browser.js("return window.queries.every(q=>q.targets===100)"),
        true
    );
    assert!(
        browser
            .js("return document.getElementById('view-title').textContent")
            .as_str()
            .unwrap()
            .contains("300 genes × 151 datasets")
    );
    for _ in 0..2 {
        browser.click("#next-rows");
        browser.wait("return !document.getElementById('download').disabled");
    }
    for _ in 0..6 {
        browser.click("#next-columns");
        browser.wait("return !document.getElementById('download').disabled");
    }
    assert!(
        browser
            .js("return document.getElementById('page-summary').textContent")
            .as_str()
            .unwrap()
            .contains("201–300")
    );
    assert_eq!(
        browser.js("return document.querySelectorAll('#heatmap td').length"),
        100
    );
    assert!(
        browser
            .js("return document.querySelector('#heatmap td').title")
            .as_str()
            .unwrap()
            .starts_with("g201 / global sum: 0")
    );
    browser.js("document.getElementById('row-position').value='1';document.getElementById('column-position').value='1'");
    browser.click("#jump-page");
    browser.wait("return !document.getElementById('download').disabled");
    assert_eq!(
        browser.js("return document.querySelectorAll('#heatmap td').length"),
        2500
    );
    browser.js("document.getElementById('row-position').value='300';document.getElementById('column-position').value='151'");
    browser.click("#jump-page");
    browser.wait("return !document.getElementById('download').disabled");
    assert_eq!(
        browser.js("return document.querySelectorAll('#heatmap td').length"),
        100
    );
    browser.js("window.exportBlob=null;const create=URL.createObjectURL;URL.createObjectURL=(file)=>{window.exportBlob=file;return create(file)}");
    browser.click("#export-selection");
    browser.wait(
        "return !document.getElementById('export-selection').disabled && window.exportBlob!==null",
    );
    let csv = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done);","args":[]}),
    );
    let csv = csv.as_str().unwrap();
    assert_eq!(csv.lines().count(), 45301);
    assert!(csv.contains("18446744073709551615"));
    assert!(csv.contains("\"300\",\"g300\",\"sample_150\",\"18446744073709551615\""));
    browser.click("#clear-export");
    browser.wait("return document.getElementById('clear-export').hidden");

    // Pure display helpers exercise exact rounding at suffix boundaries and
    // preserve full u64 counts. Changing formatting must not query more blocks.
    assert_eq!(
        browser.js("return import('./display.js').then(({formatAbundance:f})=>[f(999499n,'raw','human'),f(999500n,'raw','human'),f(1234567n,'raw','human'),f(1234567890n,'raw','human'),f(1234567890000n,'raw','human'),f(18446744073709551615n,'raw','scientific'),f(1n,'raw','scientific'),f(0n,'raw','scientific'),f(999n,'log','full')])"),
        json!(["999K","1.00M","1.23M","1.23G","1.23T","1.84e+19","1.00e+0","0.00e+0","3.00"])
    );
    let before = browser.js("return window.queries.length");
    browser.js("document.getElementById('number-format').value='scientific';document.getElementById('number-format').dispatchEvent(new Event('change'))");
    assert!(browser.js("return Array.from(document.querySelectorAll('#heatmap td')).some(cell=>cell.textContent==='1.84e+19')")==true);
    browser.js("document.getElementById('value-mode').value='log';document.getElementById('value-mode').dispatchEvent(new Event('change'))");
    assert!(browser.js("return Array.from(document.querySelectorAll('#heatmap td')).some(cell=>cell.textContent==='1.93e+1' && cell.title.includes('18446744073709551615'))")==true);
    assert_eq!(browser.js("return window.queries.length"), before);

    browser.js("document.getElementById('value-mode').value='raw';document.getElementById('hide-empty').checked=true;document.getElementById('row-sort').value='total-desc';document.getElementById('column-sort').value='total-desc'");
    browser.click("#draw");
    browser.wait("return !document.getElementById('download').disabled");
    assert!(
        browser
            .js("return document.getElementById('view-title').textContent")
            .as_str()
            .unwrap()
            .contains("226 genes × 151 datasets")
    );
    assert!(
        browser
            .js("return document.getElementById('view-details').textContent")
            .as_str()
            .unwrap()
            .contains("74 all-zero rows hidden")
    );
    assert!(
        browser
            .js("return document.querySelector('#heatmap td').title")
            .as_str()
            .unwrap()
            .starts_with("g004 / sample_026: 18446744073709551615")
    );
    assert_eq!(
        browser.js("return window.queries.slice(-25).every(q=>q.targets===100)"),
        true
    );

    // Changing name ordering reuses the summaries: only the visible page reloads.
    let before = browser.js("return window.queries.length").as_u64().unwrap();
    browser.js("document.getElementById('row-sort').value='name-asc';document.getElementById('row-sort').dispatchEvent(new Event('change'))");
    browser.click("#apply-view");
    browser.wait("return !document.getElementById('download').disabled");
    assert_eq!(
        browser.js("return window.queries.length").as_u64().unwrap() - before,
        25
    );
    assert!(
        browser
            .js("return document.querySelector('#heatmap td').title")
            .as_str()
            .unwrap()
            .starts_with("g001 / sample_026: 1")
    );

    browser.js("window.exportBlob=null");
    browser.click("#download");
    browser.wait("return window.exportBlob!==null");
    let visible = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    assert!(
        visible.as_str().unwrap().contains("18446744073709551615"),
        "scientific display must export original decoded integers"
    );
    assert!(!visible.as_str().unwrap().contains("1.84e+19"));

    // Streamed full exports follow filtering and ordering without a full matrix.
    browser.js("window.exportBlob=null");
    browser.click("#export-selection");
    browser.wait(
        "return !document.getElementById('export-selection').disabled && window.exportBlob!==null",
    );
    let filtered = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    let filtered = filtered.as_str().unwrap();
    assert_eq!(filtered.lines().count(), 226 * 151 + 1);
    assert!(filtered.contains("\"1\",\"g001\",\"sample_026\",\"1\""));
    assert!(!filtered.contains("\"g005\""));
    browser.click("#clear-export");
    browser.wait("return document.getElementById('clear-export').hidden");

    // Cancellation stops a new requested scan rather than silently publishing it.
    browser.click("#clear");
    browser.click("#select-targets");
    browser.click("#select-datasets");
    browser.click("#draw");
    browser.click("#cancel");
    browser.wait("return document.getElementById('cancel').hidden");
    assert_eq!(
        browser.js("return document.getElementById('heatmap').hidden"),
        true
    );
    assert_eq!(
        browser.js("return document.getElementById('download').disabled"),
        true
    );
}
