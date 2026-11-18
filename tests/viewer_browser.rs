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
    browser.js("window.queries=[];const send=Worker.prototype.postMessage;Worker.prototype.postMessage=function(message){window.queries.push({offset:message.offset,targets:message.request.selected.length});if(window.queryDelay){setTimeout(()=>send.call(this,message),window.queryDelay);}else return send.call(this,message);}");
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
    browser.wait("return !document.getElementById('controls').hidden && document.getElementById('controls').dataset.ready==='true'");
    assert_eq!(
        browser.js("return window.queries.length"),
        0,
        "large indexes open on a metadata overview without scanning abundance blocks"
    );
    assert_eq!(
        browser.js("return document.getElementById('main-color').value"),
        "#2563eb"
    );
    assert_eq!(
        browser.js(
            "return document.getElementById('targets-note').textContent.startsWith('300 selected')"
        ),
        true
    );
    assert_eq!(
        browser.js(
            "return document.getElementById('datasets-note').textContent.startsWith('150 selected')"
        ),
        true
    );
    assert_eq!(browser.js("return [...document.querySelectorAll('[role=tabpanel]')].filter(panel=>!panel.hidden).map(panel=>panel.id)"),json!(["overview-view"]));
    assert_eq!(
        browser.js("return document.querySelectorAll('[role=tab]').length"),
        6
    );
    assert_eq!(
        browser.js("return document.getElementById('overview-tab').getAttribute('aria-selected')"),
        "true"
    );
    // Slow work reveals progress after 500 ms; cancellation removes the bar and
    // leaves the requested update available, without publishing partial totals.
    browser.js("window.queryDelay=700");
    browser.click("#compute-overview");
    assert_eq!(
        browser.js("return document.getElementById('operation-progress').hidden"),
        true
    );
    assert_eq!(
        webdriver(
            "POST",
            &format!("{}/execute/async", browser.path),
            json!({"script":"const done=arguments[arguments.length-1];setTimeout(()=>done(document.getElementById('operation-progress').hidden),200)","args":[]}),
        ),
        true,
        "progress stays hidden during the first 200 ms"
    );
    browser.wait("return !document.getElementById('operation-progress').hidden");
    assert_eq!(
        browser.js("return document.getElementById('progress-bar').max"),
        150
    );
    browser.wait("return document.getElementById('progress-bar').value>=1");
    browser.click("#cancel");
    browser.wait("return document.getElementById('cancel').hidden");
    assert_eq!(browser.js("return document.getElementById('operation-progress').hidden && !document.getElementById('compute-overview').disabled && document.getElementById('download-overview').disabled"), true);
    browser.js("window.queryDelay=0;window.queries=[]");
    browser.click("#datasets-tab");
    browser.wait("return !document.getElementById('download-datasets').disabled && document.getElementById('cancel').hidden");
    assert_eq!(
        browser.js("return document.querySelectorAll('#datasets-table tbody tr').length"),
        25
    );
    assert_eq!(
        browser.js("return window.queries.length"),
        25,
        "dataset list scans only its visible page"
    );
    browser.js("document.getElementById('dataset-list-sort').value='total';document.getElementById('dataset-list-sort').dispatchEvent(new Event('change'))");
    assert_eq!(
        browser.js("return document.getElementById('dataset-list-sort').value"),
        "reference",
        "statistical sorting requires the complete dataset scope"
    );
    assert_eq!(
        browser.js("return window.queries.length"),
        25,
        "sorting does not trigger a hidden whole-matrix scan"
    );
    assert_eq!(
        browser.js("return window.queries.every(query=>query.targets===300)"),
        true
    );
    browser.click("#datasets-next");
    browser.wait("return document.getElementById('cancel').hidden && document.getElementById('datasets-page-label').textContent.includes('26–50')");
    assert_eq!(browser.js("return window.queries.length"), 50);
    browser.click("#datasets-table tbody button");
    browser.wait("return !document.getElementById('download-ranking').disabled && document.getElementById('cancel').hidden");
    assert_eq!(
        browser.js("return document.getElementById('ranking-tab').getAttribute('aria-selected')"),
        "true"
    );
    assert_eq!(
        browser.js("return document.getElementById('dataset-summary-title').textContent"),
        "Dataset summary · sample_026"
    );
    let before = browser.js("return window.queries.length");
    browser.click("#target-focus-tab");
    assert_eq!(
        browser.js("return window.queries.length"),
        before,
        "large target profiles wait for an explicit scan"
    );
    browser.click("#load-target-focus");
    browser.wait("return !document.getElementById('download-target-focus').disabled && document.getElementById('cancel').hidden");
    assert_eq!(
        browser.js("return window.queries.length"),
        json!(before.as_u64().unwrap() + 150)
    );
    assert_eq!(browser.js("return document.querySelector('#target-focus-cards [data-statistic=total] strong').title"),"1");
    assert_eq!(browser.js("return document.querySelector('#target-focus-cards [data-statistic=known] strong').textContent"),"150");
    assert_eq!(browser.js("return document.querySelector('#target-focus-table [data-statistic=value]').textContent"),"0");
    let before = browser.js("return window.queries.length");
    browser.click("#target-focus-next");
    browser.wait("return document.getElementById('cancel').hidden && document.getElementById('target-focus-page-label').textContent.includes('26–50')");
    assert_eq!(
        browser.js("return window.queries.length"),
        json!(before.as_u64().unwrap() + 25)
    );
    assert_eq!(
        browser.js(
            "return document.querySelector('#target-focus-table [data-statistic=fold]').textContent"
        ),
        "150"
    );
    assert_eq!(
        browser.js("return document.querySelectorAll('#target-focus-table tbody tr').length"),
        25
    );
    let before = browser.js("return window.queries.length");
    browser.js("document.getElementById('color-preset').value='#7c3aed';document.getElementById('color-preset').dispatchEvent(new Event('change'))");
    assert_eq!(
        browser.js(
            "return document.querySelector('#ranking-chart rect[data-share]').getAttribute('fill')"
        ),
        "#7c3aed"
    );
    assert_eq!(
        browser.js("return document.documentElement.style.getPropertyValue('--accent')"),
        "#7c3aed"
    );
    assert_eq!(
        browser.js("return window.queries.length"),
        before,
        "color changes reuse all cached data"
    );
    browser.js("document.getElementById('main-color').value='#ffff00';document.getElementById('main-color').dispatchEvent(new Event('input'))");
    assert_eq!(
        browser.js("return document.documentElement.style.getPropertyValue('--on-accent')"),
        "#172b4d"
    );
    assert_eq!(
        browser.js("return document.getElementById('color-preset').value"),
        "custom"
    );
    assert_eq!(browser.js("return window.queries.length"), before);
    browser.js("document.getElementById('color-preset').value='#2563eb';document.getElementById('color-preset').dispatchEvent(new Event('change'));window.queries=[]");
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

    browser.js("document.getElementById('value-mode').value='raw';document.getElementById('hide-empty').checked=true;document.getElementById('row-sort').value='total-desc';document.getElementById('column-sort').value='total-desc';document.getElementById('hide-empty').dispatchEvent(new Event('change'))");
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
            .contains("74 empty rows hidden")
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

    // Legacy indexes remain usable in raw mode; unavailable totals fail before any scan.
    let before = browser.js("return window.queries.length");
    browser.js("document.getElementById('normalization').value='cpm';document.getElementById('normalization').dispatchEvent(new Event('change'))");
    browser.click("#draw");
    browser.wait("return document.getElementById('status').textContent.includes('Dataset k-mer totals are unavailable')");
    assert_eq!(browser.js("return window.queries.length"), before);
    browser.js("document.getElementById('normalization').value='none';document.getElementById('normalization').dispatchEvent(new Event('change'))");

    // The reference must include the dataset on page two where g001 alone is expressed.
    browser.js("document.getElementById('value-mode').value='fold';document.getElementById('value-mode').dispatchEvent(new Event('change'));document.getElementById('row-sort').value='reference';document.getElementById('column-sort').value='reference';document.getElementById('hide-empty').checked=false;document.getElementById('number-format').value='full';document.getElementById('min').value='';document.getElementById('max').value=''");
    browser.click("#draw");
    browser.wait("return !document.getElementById('download').disabled");
    assert_eq!(
        browser.js("return document.querySelector('#heatmap td').textContent"),
        "0"
    );
    assert_eq!(browser.js("return document.querySelector('#heatmap td').title.includes('150 known selected real datasets')"),true);
    browser.js("document.getElementById('column-position').value='26'");
    browser.click("#jump-page");
    browser.wait("return !document.getElementById('download').disabled");
    assert_eq!(
        browser.js("return document.querySelector('#heatmap td').textContent"),
        "150"
    );
    assert_eq!(
        browser.js("return document.querySelectorAll('#heatmap td').length"),
        2500
    );
    browser.js("document.getElementById('value-mode').value='raw';document.getElementById('value-mode').dispatchEvent(new Event('change'))");

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

    let toy = root.join("examples/viewer/toy.eai");
    webdriver(
        "POST",
        &format!("{}/element/{element}/value", browser.path),
        json!({"text":toy.to_str().unwrap()}),
    );
    browser.wait("return document.getElementById('target-count').textContent==='6' && document.getElementById('controls').dataset.ready==='true'");
    assert_eq!(
        browser.js("return document.getElementById('overview-tab').getAttribute('aria-selected')"),
        "true"
    );
    assert_eq!(
        browser.js(
            "return document.querySelector('#overview-cards [data-statistic=total] strong').title"
        ),
        "19278"
    );
    assert_eq!(
        browser.js(
            "return document.querySelector('#overview-cards [data-statistic=known] strong').title"
        ),
        "24"
    );
    assert_eq!(
        browser.js(
            "return document.querySelector('#overview-cards [data-statistic=nonzero] strong').title"
        ),
        "16"
    );
    assert_eq!(browser.js("return document.querySelector('#overview-cards [data-statistic=median] strong').textContent"),"85");
    assert_eq!(browser.js("return document.getElementById('compute-overview').disabled && !document.getElementById('compute-overview').classList.contains('needs-update')"),true);
    let before = browser.js("return window.queries.length");
    browser.click("#compute-overview");
    assert_eq!(
        browser.js("return window.queries.length"),
        before,
        "current summaries cannot be recomputed with an inactive button"
    );
    browser.js("document.getElementById('overview-threshold').value='100';document.getElementById('overview-threshold').dispatchEvent(new Event('input',{bubbles:true}))");
    assert_eq!(browser.js("return !document.getElementById('compute-overview').disabled && document.getElementById('compute-overview').classList.contains('needs-update')"),true);
    browser.js("document.getElementById('overview-threshold').value='0';document.getElementById('overview-threshold').dispatchEvent(new Event('input',{bubbles:true}))");
    assert_eq!(
        browser.js("return document.getElementById('compute-overview').disabled"),
        true,
        "reverting a pending parameter restores the current state"
    );
    assert_eq!(browser.js("return !document.body.textContent.includes('private by design')&&!document.body.textContent.includes('FROM SEQUENCE TO SIGNAL')&&!document.body.textContent.includes('Your counts, in context')&&!document.body.textContent.includes('LOCAL DECODING')"),true);
    browser.click("#datasets-tab");
    browser.wait("return !document.getElementById('download-datasets').disabled && document.getElementById('cancel').hidden");
    assert_eq!(browser.js("return document.getElementById('compute-dataset-page').disabled&&document.getElementById('compute-datasets-all').disabled"),true);
    assert_eq!(browser.js("return document.querySelector('#datasets-table tbody [data-statistic=median]').textContent"),"51");
    browser.js("document.getElementById('dataset-list-sort').value='total';document.getElementById('dataset-list-sort').dispatchEvent(new Event('change'))");
    assert_eq!(
        browser.js("return document.querySelector('#datasets-table tbody button').textContent"),
        "sample_B2"
    );
    browser.js("document.getElementById('dataset-list-sort').value='reference';document.getElementById('dataset-list-sort').dispatchEvent(new Event('change'))");
    browser.click("#datasets-table tbody button");
    browser.wait("return !document.getElementById('download-ranking').disabled && document.getElementById('cancel').hidden");
    assert_eq!(browser.js("return document.querySelector('#dataset-summary-cards [data-statistic=total] strong').title"),"4862");
    assert_eq!(browser.js("return document.querySelector('#dataset-summary-cards [data-statistic=median] strong').textContent"),"51");
    assert_eq!(browser.js("return document.querySelector('#dataset-summary-cards [data-statistic=nonzero] strong').textContent"),"4");
    assert_eq!(
        browser.js("return document.getElementById('draw-ranking').disabled"),
        true
    );
    assert_eq!(browser.js("return document.getElementById('targets-note').textContent.startsWith('6 selected') && document.getElementById('datasets-note').textContent.startsWith('4 selected')"),true);
    browser.click("#target-summary-tab");
    browser.wait("return document.getElementById('cancel').hidden");
    assert_eq!(
        browser.js("return document.getElementById('target-summary-panel').hidden"),
        false
    );
    browser.js("window.exportBlob=null");
    browser.click("#download-dataset-summary");
    browser.wait("return window.exportBlob!==null");
    let dataset_csv = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    let mut reader = csv::Reader::from_reader(dataset_csv.as_str().unwrap().as_bytes());
    let header = reader.headers().unwrap().clone();
    let median = header.iter().position(|h| h == "median").unwrap();
    let row = reader.records().next().unwrap().unwrap();
    assert_eq!(&row[median], "51");
    browser.click("#summary-table tbody button");
    browser.wait("return !document.getElementById('download-target-focus').disabled && document.getElementById('cancel').hidden");
    assert_eq!(
        browser
            .js("return document.getElementById('target-focus-tab').getAttribute('aria-selected')"),
        "true"
    );
    assert_eq!(browser.js("return document.querySelector('#target-focus-cards [data-statistic=total] strong').title"),"7310");
    assert_eq!(browser.js("return document.querySelector('#target-focus-cards [data-statistic=median] strong').textContent"),"1,581");
    assert_eq!(
        browser.js("return document.getElementById('load-target-focus').disabled"),
        true
    );
    assert_eq!(
        browser.js(
            "return document.querySelector('#target-focus-table [data-statistic=fold]').textContent"
        ),
        "2.23256"
    );
    browser.js("window.exportBlob=null");
    browser.click("#download-target-page");
    browser.wait("return window.exportBlob!==null");
    let profile_csv = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    let mut reader = csv::Reader::from_reader(profile_csv.as_str().unwrap().as_bytes());
    assert_eq!(reader.records().count(), 4);
    browser.js("document.getElementById('target-list').value='#1'");
    browser.click("#apply-targets");
    assert_eq!(
        browser.js(
            "return document.getElementById('targets-note').textContent.startsWith('1 selected')"
        ),
        true,
        "pasting an ID replaces the default all-target selection"
    );
    assert_eq!(
        browser.js("return document.getElementById('dataset-list-sort').value"),
        "reference"
    );
    assert_eq!(
        browser.js(
            "return document.querySelector('#overview-cards [data-statistic=total] strong').title"
        ),
        "19278",
        "global overview stays independent of target selection"
    );
    browser.click("#compute-summary");
    browser.wait("return !document.getElementById('download-summary').disabled");
    assert_eq!(
        browser.js("return document.querySelectorAll('#summary-table tbody tr').length"),
        1
    );
    assert_eq!(
        browser.js("return document.querySelector('#summary-table tbody th').textContent"),
        "GENE_A"
    );
    assert_eq!(browser.js("return document.querySelector('#summary-table tbody [data-statistic=total]').textContent"),"7,310");
    assert_eq!(
        browser.js("return document.getElementById('compute-summary').disabled"),
        true
    );
    // Fold-change reference spans the complete dataset selection, excludes aggregate
    // vectors, includes known zeros, and stays independent of bounds and paging.
    browser.js("document.getElementById('targets-search').value='';document.getElementById('datasets-search').value='';document.getElementById('normalization').value='none';document.getElementById('normalization').dispatchEvent(new Event('change'));document.getElementById('value-mode').value='fold';document.getElementById('value-mode').dispatchEvent(new Event('change'));document.getElementById('hide-empty').checked=false;document.getElementById('row-sort').value='reference';document.getElementById('column-sort').value='reference';document.getElementById('min').value='';document.getElementById('max').value='';document.getElementById('number-format').value='full'");
    browser.click("#select-targets");
    browser.click("#select-datasets");
    browser.click("#draw");
    browser.wait("return !document.getElementById('download').disabled");
    assert_eq!(
        browser.js("return document.querySelector('#heatmap td').textContent"),
        "2.23256"
    );
    assert_eq!(browser.js("return document.querySelector('#heatmap td').title.includes('target mean: 1827.5 (4 known')"),true);
    assert_eq!(browser.js("return document.querySelectorAll('#heatmap tbody tr')[4].querySelector('td').textContent"),"1");
    assert_eq!(browser.js("return document.querySelectorAll('#heatmap tbody tr')[5].querySelector('td').textContent"),"—");
    assert_eq!(browser.js("return document.getElementById('draw').disabled && document.getElementById('apply-view').disabled"),true);
    assert_eq!(browser.js("return [...document.querySelector('#heatmap tbody tr').querySelectorAll('td')].slice(4).every(td=>td.textContent==='—')"),true);
    let before = browser.js("return window.queries.length");
    browser.js("document.getElementById('value-mode').value='log-fold';document.getElementById('value-mode').dispatchEvent(new Event('change'))");
    assert_eq!(
        browser.js("return window.queries.length"),
        before,
        "FC/log FC switches reuse reference means"
    );
    assert_eq!(
        browser.js("return document.querySelector('#heatmap td').textContent"),
        "1.1587"
    );
    assert_eq!(browser.js("return document.querySelectorAll('#heatmap tbody tr')[2].querySelector('td').textContent"),"−∞");
    assert_eq!(browser.js("return document.querySelectorAll('#heatmap tbody tr')[4].querySelector('td').textContent"),"0");
    assert_eq!(
        browser.js("return document.getElementById('scale').disabled"),
        true
    );
    browser.js("document.getElementById('min').value='100';document.getElementById('min').dispatchEvent(new Event('change'))");
    assert_eq!(browser.js("return window.queries.length"), before);
    assert_eq!(
        browser.js("return document.querySelector('#heatmap td').textContent"),
        "1.1587",
        "bounds must not alter means"
    );
    browser.js("document.getElementById('min').value='';document.getElementById('min').dispatchEvent(new Event('change'));window.exportBlob=null");
    browser.click("#download");
    browser.wait("return window.exportBlob!==null");
    let fold_csv = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    let mut reader = csv::Reader::from_reader(fold_csv.as_str().unwrap().as_bytes());
    let rows: Vec<_> = reader.records().map(Result::unwrap).collect();
    assert!((rows[0][2].parse::<f64>().unwrap() - (4080.0_f64 / 1827.5).log2()).abs() < 1e-12);
    assert_eq!(&rows[2][2], "-Infinity");
    assert_eq!(&rows[5][2], "");
    assert_eq!(&rows[0][6], "");
    browser.js("window.exportBlob=null");
    browser.click("#export-selection");
    browser.wait(
        "return !document.getElementById('export-selection').disabled && window.exportBlob!==null",
    );
    let fold_csv = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    let mut reader = csv::Reader::from_reader(fold_csv.as_str().unwrap().as_bytes());
    assert_eq!(reader.headers().unwrap().get(3), Some("log2_fold_change"));
    let rows: Vec<_> = reader.records().map(Result::unwrap).collect();
    assert_eq!(rows.len(), 20);
    assert!(
        rows.iter()
            .all(|row| row[2].starts_with("sample_") && &row[1] != "GENE_EMPTY")
    );
    browser.click("#clear-export");
    browser.wait("return document.getElementById('clear-export').hidden");
    // Changed dataset selection recomputes the mean, after normalization.
    browser.click("#clear");
    browser.js("document.getElementById('target-list').value='GENE_A';document.getElementById('dataset-list').value=['sample_A1','sample_A2'].join(String.fromCharCode(10));document.getElementById('normalization').value='cpm';document.getElementById('normalization').dispatchEvent(new Event('change'));document.getElementById('value-mode').value='fold'");
    browser.click("#apply-targets");
    browser.click("#apply-datasets");
    browser.click("#draw");
    browser.wait("return !document.getElementById('download').disabled");
    browser.js("window.exportBlob=null");
    browser.click("#download");
    browser.wait("return window.exportBlob!==null");
    let fold_csv = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    let mut reader = csv::Reader::from_reader(fold_csv.as_str().unwrap().as_bytes());
    let row = reader.records().next().unwrap().unwrap();
    let expected = (4080.0 / 4862.0) / ((4080.0 / 4862.0 + 3060.0 / 3842.0) / 2.0);
    assert!((row[2].parse::<f64>().unwrap() - expected).abs() < 1e-12);
    assert_eq!(browser.js("return import('./display.js').then(({foldChange:f,formatAbundance:p})=>[f(0n,10n,2),f(0n,10n,2,'log-fold')===-Infinity,f(0n,0n,2),f(null,10,1),f(3.5,3.5,1,'log-fold'),p(-Infinity),f(18446744073709551615n,36893488147419103230n,2)])"),json!([0,true,null,null,0,"−∞",1]));
    browser.click("#clear");
    browser.js("document.getElementById('target-list').value='GENE_A';document.getElementById('dataset-list').value='sample_A1';document.getElementById('hide-empty').checked=false;document.getElementById('row-sort').value='reference';document.getElementById('column-sort').value='reference';document.getElementById('min').value='';document.getElementById('max').value='';document.getElementById('number-format').value='scientific';document.getElementById('value-mode').value='raw'");
    browser.click("#apply-targets");
    browser.click("#apply-datasets");
    browser.js("document.getElementById('normalization').value='length-cpm';document.getElementById('normalization').dispatchEvent(new Event('change'))");
    let before = browser.js("return window.queries.length").as_u64().unwrap();
    browser.click("#draw");
    browser.wait("return !document.getElementById('download').disabled");
    assert_eq!(
        browser.js("return window.queries.length").as_u64().unwrap() - before,
        1,
        "normalization loads only the requested dataset and target"
    );
    assert_eq!(
        browser.js("return document.querySelector('#heatmap td').textContent"),
        "2.47e+4"
    );
    assert!(
        browser
            .js("return document.querySelector('#heatmap td').title.includes('raw count: 4080')")
            == true
    );
    browser.js("window.exportBlob=null");
    browser.click("#download");
    browser.wait("return window.exportBlob!==null");
    let normalized = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    let mut reader = csv::Reader::from_reader(normalized.as_str().unwrap().as_bytes());
    let row = reader.records().next().unwrap().unwrap();
    assert!(
        (row[2].parse::<f64>().unwrap() - 4080.0 / 34.0 / 4862.0 * 1e6).abs() < 1e-8,
        "full dataset denominator is preserved despite selecting one target"
    );
    browser.js("window.exportBlob=null");
    browser.click("#export-selection");
    browser.wait(
        "return !document.getElementById('export-selection').disabled && window.exportBlob!==null",
    );
    let normalized = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    let mut reader = csv::Reader::from_reader(normalized.as_str().unwrap().as_bytes());
    let row = reader.records().next().unwrap().unwrap();
    assert!((row[3].parse::<f64>().unwrap() - 4080.0 / 34.0 / 4862.0 * 1e6).abs() < 1e-8);
    browser.click("#clear-export");
    browser.wait("return document.getElementById('clear-export').hidden");

    // Pure math covers fractional/zero totals and target lengths, and normalized display.
    assert_eq!(browser.js("return import('./display.js').then(({normalize:n,formatAbundance:f})=>{const d={kind:'dataset',total_kmers:{numerator:'4000001',scale:1000000}};return [n(1n,{kmer_length:10},d,'coverage'),Math.round(n(1n,{kmer_length:10},d,'length-cpm')*1e8)/1e8,n(0n,{kmer_length:0},{kind:'dataset',total_kmers:{numerator:'0',scale:1}},'length-cpm'),f(0.125,'raw','full'),f(0.125,'raw','scientific')];})"), json!([0.1,24999.99375,0,"0.125","1.25e-1"]));

    // Sorting totals uses normalized values, not raw counts; statistics are rejected.
    browser.click("#clear");
    browser.click("#select-targets");
    browser.js("document.getElementById('datasets-search').value='sample_';document.getElementById('datasets-search').dispatchEvent(new Event('input'));document.getElementById('normalization').value='cpm';document.getElementById('normalization').dispatchEvent(new Event('change'));document.getElementById('row-sort').value='total-desc';document.getElementById('column-sort').value='total-desc'");
    browser.click("#select-datasets");
    browser.click("#draw");
    browser.wait("return !document.getElementById('download').disabled");
    assert_eq!(
        browser.js("return document.querySelector('#heatmap tbody th').textContent"),
        "GENE_A",
        "normalization reverses raw gene B > A total ordering"
    );
    browser.js("document.getElementById('dataset-list').value='mean'");
    browser.click("#apply-datasets");
    browser.click("#draw");
    browser.wait("return document.getElementById('status').textContent.includes('Remove raw statistic vectors')");

    // On-demand summaries exclude global/statistic entries, include zeros, and ignore heatmap bounds.
    browser.js("document.getElementById('summary-threshold').value='100';document.getElementById('number-format').value='full';document.getElementById('min').value='999999';document.getElementById('summary-percentiles').value='0,10,50,90,100'");
    let before = browser.js("return window.queries.length").as_u64().unwrap();
    browser.click("#compute-summary");
    browser.wait("return !document.getElementById('download-summary').disabled");
    assert_eq!(
        browser.js("return window.queries.length").as_u64().unwrap() - before,
        4,
        "all-dataset summaries must scan the four real datasets once, excluding global/statistics"
    );
    assert_eq!(browser.js("return [...document.querySelectorAll('#summary-table tbody tr')].find(row=>row.querySelector('th').textContent==='GENE_A').querySelector('[data-statistic=total]').textContent"),"7,310");
    assert_eq!(browser.js("return [...document.querySelectorAll('#summary-table tbody tr')].find(row=>row.querySelector('th').textContent==='GENE_A').querySelector('[data-statistic=above]').textContent"),"3");
    assert_eq!(browser.js("return [...document.querySelectorAll('#summary-table tbody tr')].find(row=>row.querySelector('th').textContent==='GENE_A').querySelector('[data-statistic=mean]').textContent"),"1,827.5");
    assert_eq!(browser.js("return [...document.querySelectorAll('#summary-table tbody tr')].find(row=>row.querySelector('th').textContent==='GENE_A').querySelector('[data-statistic=median]').textContent"),"1,581");
    assert_eq!(browser.js("return [...document.querySelectorAll('#summary-table tbody tr')].find(row=>row.querySelector('th').textContent==='GENE_A').querySelector('[data-statistic=p10]').textContent"),"78.2");
    let before = browser.js("return window.queries.length");
    browser.js("document.getElementById('summary-sort').value='total';document.getElementById('summary-sort').dispatchEvent(new Event('change'))");
    assert_eq!(
        browser.js("return document.querySelector('#summary-table tbody th').textContent"),
        "GENE_B"
    );
    assert_eq!(
        browser.js("return window.queries.length"),
        before,
        "sorting summary rows uses cached summaries"
    );
    assert_eq!(browser.js("return document.querySelector('#summary-table tbody [data-statistic=total]').style.backgroundColor!==''"),true);
    browser.js("document.getElementById('summary-color').checked=false;document.getElementById('summary-color').dispatchEvent(new Event('change'));document.getElementById('summary-sort').value='name';document.getElementById('summary-sort').dispatchEvent(new Event('change'))");
    assert_eq!(
        browser.js("return document.querySelector('#summary-table tbody th').textContent"),
        "HOUSEKEEPING"
    );
    assert_eq!(browser.js("return document.querySelector('#summary-table tbody [data-statistic=total]').style.backgroundColor"),"");
    browser.js("document.getElementById('summary-direction').value='asc';document.getElementById('summary-direction').dispatchEvent(new Event('change'))");
    assert_eq!(
        browser.js("return document.querySelector('#summary-table tbody th').textContent"),
        "GENE_A"
    );
    browser.js("document.getElementById('summary-direction').value='desc';document.getElementById('summary-sort').value='total';document.getElementById('summary-sort').dispatchEvent(new Event('change'));document.getElementById('summary-color').checked=true;document.getElementById('summary-color').dispatchEvent(new Event('change'))");
    assert_eq!(
        browser.js("return window.queries.length"),
        before,
        "presentation and name sorting do not read abundance blocks"
    );
    browser.js("window.exportBlob=null");
    browser.click("#download-summary");
    browser.wait("return window.exportBlob!==null");
    let summary_csv = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    let mut reader = csv::Reader::from_reader(summary_csv.as_str().unwrap().as_bytes());
    let header = reader.headers().unwrap().clone();
    let total = header.iter().position(|h| h == "total").unwrap();
    let mean = header.iter().position(|h| h == "mean").unwrap();
    let records = reader.records().map(Result::unwrap).collect::<Vec<_>>();
    assert_eq!(records.len(), 6);
    assert_eq!(&records[0][1], "GENE_B");
    assert_eq!(&records[0][total], "8330");
    assert_eq!(
        &records.iter().find(|r| &r[1] == "GENE_A").unwrap()[mean],
        "1827.5"
    );

    // A selected-only summary uses the chosen normalization and full library denominator.
    browser.click("#clear");
    browser.click("#select-targets");
    browser.js("document.getElementById('dataset-list').value='sample_A1';document.getElementById('summary-scope').value='selected';document.getElementById('summary-values').value='view'");
    browser.click("#apply-datasets");
    browser.click("#compute-summary");
    browser.wait("return !document.getElementById('download-summary').disabled");
    assert_eq!(browser.js("return [...document.querySelectorAll('#summary-table tbody tr')].find(row=>row.querySelector('th').textContent==='GENE_A').querySelector('[data-statistic=known]').textContent"),"1");
    browser.js("window.exportBlob=null");
    browser.click("#download-summary");
    browser.wait("return window.exportBlob!==null");
    let normalized_summary = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    let mut reader = csv::Reader::from_reader(normalized_summary.as_str().unwrap().as_bytes());
    let header = reader.headers().unwrap().clone();
    let total = header.iter().position(|h| h == "total").unwrap();
    let normalized_gene = reader
        .records()
        .map(Result::unwrap)
        .find(|row| &row[1] == "GENE_A")
        .unwrap();
    assert!((normalized_gene[total].parse::<f64>().unwrap() - 4080.0 / 4862.0 * 1e6).abs() < 1e-8);
    browser.js("document.getElementById('summary-scope').value='all';document.getElementById('summary-values').value='raw';document.getElementById('summary-percentiles').value='101';document.getElementById('summary-percentiles').dispatchEvent(new Event('change'))");
    let before = browser.js("return window.queries.length");
    browser.click("#compute-summary");
    browser.wait("return document.getElementById('status').textContent.includes('Percentiles must be between')");
    assert_eq!(
        browser.js("return window.queries.length"),
        before,
        "invalid settings must fail before scanning"
    );
    browser.js("document.getElementById('summary-percentiles').value='0,10,50,90,100'");

    // Zero exclusion changes distribution statistics, while observed total remains zero.
    browser.js("document.getElementById('summary-exclude-zeros').checked=true;document.getElementById('summary-exclude-zeros').dispatchEvent(new Event('change'))");
    browser.click("#compute-summary");
    browser.wait("return !document.getElementById('download-summary').disabled");
    assert_eq!(browser.js("return [...document.querySelectorAll('#summary-table tbody tr')].find(row=>row.querySelector('th').textContent==='GENE_EMPTY').querySelector('[data-statistic=mean]').textContent"),"—");

    // Independent summary math: large integer precision, strict threshold, percentile definitions,
    // missing values, and deterministic bounded sampling of a large distribution.
    assert_eq!(browser.js("return import('./statistics.js').then(({Distribution:D,compareStatistics:c})=>{const d=new D({threshold:18446744073709551614n});d.add(18446744073709551613n);d.add(18446744073709551615n);d.add(null);const r=d.result([25,50]);return [r.total.toString(),r.mean,r.median,r.percentiles[25],r.above,r.known,r.missing,r.stddev,c(r.mean,18446744073709551613n)];})"),json!(["36893488147419103228","18446744073709551614","18446744073709551614","18446744073709551613.5",1,2,1,1,1]));
    assert_eq!(browser.js("return import('./statistics.js').then(({Distribution:D})=>{const a=new D({method:'sample',sampleSize:256});const b=new D({method:'sample',sampleSize:256});for(let i=0;i<100000;i++){a.add(BigInt(i));b.add(BigInt(i));}const x=a.result([0,50,100]);const y=b.result([50]);return [a.sample.length,a.histogram.size,x.total.toString(),x.quantiles_approximate,x.percentiles[0].toString(),x.percentiles[100].toString(),x.percentiles[50]===y.percentiles[50]];})"),json!([256,0,"4999950000",true,"0","99999",true]));
    assert_eq!(browser.js("return import('./statistics.js').then(({Distribution:D})=>{const d=new D();[0n,10n,20n,30n].forEach(x=>d.add(x));const rank=d.result([25,50],'nearest-rank');const linear=d.result([25,50],'linear');return [rank.percentiles[25].toString(),rank.percentiles[50].toString(),rank.median,linear.percentiles[25],linear.percentiles[50]];})"),json!(["0","10","15","7.5","15"]));

    // Rankings scan one dataset, sum every target, and retain only a bounded top N.
    let before = browser.js("return window.queries.length");
    browser.click("#ranking-tab");
    browser.wait("return !document.getElementById('download-ranking').disabled && document.getElementById('cancel').hidden");
    assert_eq!(
        browser.js("return window.queries.length"),
        json!(before.as_u64().unwrap() + 1)
    );
    let before = browser.js("return window.queries.length");
    browser.js("document.getElementById('ranking-limit').value='2';document.getElementById('ranking-format').value='full';document.getElementById('ranking-other').checked=false;document.getElementById('ranking-limit').dispatchEvent(new Event('change'))");
    browser.click("#draw-ranking");
    browser.wait("return !document.getElementById('download-ranking').disabled");
    assert_eq!(
        browser.js("return window.queries.length"),
        json!(before.as_u64().unwrap() + 1)
    );
    assert_eq!(
        browser.js(
            "return [...document.querySelectorAll('#ranking-chart g')].map(g=>g.dataset.target)"
        ),
        json!(["0", "4"])
    );
    let share = browser.js("return Number(document.querySelector('#ranking-chart rect[data-share]').dataset.share)").as_f64().unwrap();
    assert!((share - 4080.0 / 4862.0 * 100.0).abs() < 1e-10);
    assert_eq!(
        browser
            .js("return document.getElementById('ranking-details').textContent.includes('97.90%')"),
        true
    );
    let before = browser.js("return window.queries.length");
    browser.js("document.getElementById('ranking-other').checked=true;document.getElementById('ranking-other').dispatchEvent(new Event('change'));document.getElementById('ranking-axis').value='fixed';document.getElementById('ranking-axis').dispatchEvent(new Event('change'))");
    assert_eq!(browser.js("return window.queries.length"), before);
    assert_eq!(browser.js("return document.querySelector('#ranking-chart g[data-target=other]').textContent.includes('102')"),true);
    browser.js("window.exportBlob=null");
    browser.click("#download-ranking");
    browser.wait("return window.exportBlob!==null");
    let ranking_csv = webdriver(
        "POST",
        &format!("{}/execute/async", browser.path),
        json!({"script":"const done=arguments[arguments.length-1];window.exportBlob.text().then(done)","args":[]}),
    );
    let mut reader = csv::Reader::from_reader(ranking_csv.as_str().unwrap().as_bytes());
    let header = reader.headers().unwrap().clone();
    let denominator = header.iter().position(|h| h == "denominator").unwrap();
    let value = header.iter().position(|h| h == "value").unwrap();
    let rows: Vec<_> = reader.records().map(Result::unwrap).collect();
    assert_eq!(rows.len(), 3);
    assert_eq!(&rows[0][denominator], "4862");
    assert_eq!(&rows[0][value], "4080");
    assert_eq!(&rows[2][value], "102");
    browser.js("window.exportBlob=null");
    browser.click("#download-ranking-svg");
    browser.wait("return window.exportBlob!==null");
    assert_eq!(browser.js("return window.exportBlob.text().then(s=>s.includes('http://www.w3.org/2000/svg')&&s.includes('GENE_A')&&s.includes('83.92%'))"),true);
    browser.click("#ranking-next");
    browser.wait("return !document.getElementById('download-ranking').disabled");
    assert_eq!(
        browser.js("return document.getElementById('ranking-title').textContent"),
        "Most expressed genes · sample_A2"
    );
    assert_eq!(
        browser
            .js("return document.getElementById('ranking-details').textContent.includes('3,842')"),
        true
    );

    // A candidate subset must not inflate percentages by shrinking the total.
    browser.click("#clear");
    browser.js("document.getElementById('target-list').value='GENE_B'");
    browser.click("#apply-targets");
    browser.js("document.getElementById('ranking-targets').value='selected';document.getElementById('ranking-targets').dispatchEvent(new Event('change'));document.getElementById('ranking-other').checked=false");
    browser.click("#draw-ranking");
    browser.wait("return !document.getElementById('download-ranking').disabled");
    assert_eq!(
        browser.js(
            "return [...document.querySelectorAll('#ranking-chart g')].map(g=>g.dataset.target)"
        ),
        json!(["1"])
    );
    let share=browser.js("return Number(document.querySelector('#ranking-chart rect[data-share]').dataset.share)").as_f64().unwrap();
    assert!((share - 102.0 / 3842.0 * 100.0).abs() < 1e-10);
    browser.js("document.getElementById('ranking-denominator').value='input';document.getElementById('ranking-denominator').dispatchEvent(new Event('change'))");
    browser.click("#draw-ranking");
    browser.wait("return !document.getElementById('download-ranking').disabled");
    assert_eq!(browser.js("return document.getElementById('ranking-details').textContent.includes('full dataset k-mer total')"),true);
    browser.js("document.getElementById('ranking-values').value='cpm';document.getElementById('ranking-values').dispatchEvent(new Event('change'))");
    assert_eq!(
        browser.js(
            "return document.querySelector('#ranking-denominator option[value=input]').disabled"
        ),
        true
    );
    browser.click("#draw-ranking");
    browser.wait("return !document.getElementById('download-ranking').disabled");
    let share=browser.js("return Number(document.querySelector('#ranking-chart rect[data-share]').dataset.share)").as_f64().unwrap();
    assert!((share - 102.0 / 3842.0 * 100.0).abs() < 1e-10);
    browser.js("document.getElementById('ranking-name').value='sample_B1'");
    browser.click("#ranking-find");
    browser.wait("return !document.getElementById('download-ranking').disabled");
    assert_eq!(
        browser.js("return document.getElementById('ranking-title').textContent"),
        "Most expressed genes · sample_B1"
    );
    browser.js("document.getElementById('ranking-limit').value='0';document.getElementById('ranking-limit').dispatchEvent(new Event('input',{bubbles:true}))");
    let before = browser.js("return window.queries.length");
    browser.click("#draw-ranking");
    browser.wait("return document.getElementById('status').textContent.includes('Top N must')");
    assert_eq!(browser.js("return window.queries.length"), before);
    browser.js("document.getElementById('ranking-limit').value='2'");
    browser.click("#clear");
    browser.js("document.getElementById('target-list').value='GENE_EMPTY'");
    browser.click("#apply-targets");
    browser.click("#draw-ranking");
    browser.wait("return !document.getElementById('download-ranking').disabled");
    assert_eq!(browser.js("return document.getElementById('ranking-chart').textContent.includes('No positive expression')"),true);
    browser.js("document.getElementById('ranking-datasets').value='selected';document.getElementById('ranking-datasets').dispatchEvent(new Event('change'))");
    let before = browser.js("return window.queries.length");
    browser.click("#draw-ranking");
    browser.wait("return document.getElementById('status').textContent.includes('Select at least one real dataset')");
    assert_eq!(browser.js("return window.queries.length"), before);
    browser.js("document.getElementById('ranking-datasets').value='all';document.getElementById('ranking-datasets').dispatchEvent(new Event('change'));document.getElementById('ranking-targets').value='all'");

    // Independent heap check against a full sort; integer totals must stay exact.
    assert_eq!(browser.js("return import('./ranking.js').then(({Ranking:R,share})=>{const r=new R(17),values=[];for(let id=0;id<100000;id++){const value=18446744073709551615n-BigInt((id*7919)%10007);values.push({id,value});r.add(id,value);}const expected=values.sort((a,b)=>a.value>b.value?-1:a.value<b.value?1:a.id-b.id).slice(0,17);const result=r.result();return [r.heap.length,result.rows.every((row,i)=>row.id===expected[i].id),result.total.toString(),share(0n,0n)];})"),json!([17,true,"1844674407370954661195082",null]));
    browser.click("#comparison-tab");

    // Transcript packing and the same statistics UI preserve fractional and missing observations.
    let transcripts = dir.join("transcript-results");
    fs::create_dir(&transcripts).unwrap();
    let header = "transcript_id,transcript_name,gene_name,observable_kmers,abundance,assigned_observations\n";
    fs::write(
        transcripts.join("A.csv"),
        format!("{header}1,T1,G,10,1.25,12.5\n2,T2,G,20,,\n3,T3,H,5,0,0\n"),
    )
    .unwrap();
    fs::write(
        transcripts.join("B.csv"),
        format!("{header}1,T1,G,10,2.5,25\n2,T2,G,20,,\n3,T3,H,5,3,15\n"),
    )
    .unwrap();
    fs::write(
        transcripts.join("global.csv"),
        format!("{header}1,T1,G,10,3.75,37.5\n2,T2,G,20,,\n3,T3,H,5,3,15\n"),
    )
    .unwrap();
    fs::write(transcripts.join("manifest.json"),serde_json::to_vec(&json!({"version":1,"format":"isoform-abundance","datasets":[{"name":"A","output":"A.csv"},{"name":"B","output":"B.csv"}],"global":"global.csv","source_coverage":{"complete":true}})).unwrap()).unwrap();
    let transcript_index = dir.join("transcripts.eai");
    let result = Command::new(env!("CARGO_BIN_EXE_expresso"))
        .args(["pack", "--input"])
        .arg(&transcripts)
        .arg("--output")
        .arg(&transcript_index)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let before = browser.js("return window.queries.length");
    webdriver(
        "POST",
        &format!("{}/element/{element}/value", browser.path),
        json!({"text":transcript_index.to_str().unwrap()}),
    );
    browser.wait("return document.getElementById('target-kind').textContent==='TRANSCRIPTS' && document.getElementById('controls').dataset.ready==='true'");
    assert_eq!(
        browser.js("return window.queries.length"),
        json!(before.as_u64().unwrap() + 2),
        "opening a tiny transcript matrix computes its whole-matrix overview"
    );
    assert_eq!(browser.js("return document.querySelector('#overview-cards [data-statistic=missing] strong').textContent"),"2");
    browser.click("#select-targets");
    browser.js("document.getElementById('summary-threshold').value='10';document.getElementById('summary-exclude-zeros').checked=false;document.getElementById('min').value='';document.getElementById('summary-sort').value='reference'");
    browser.click("#compute-summary");
    browser.wait("return !document.getElementById('download-summary').disabled");
    assert_eq!(browser.js("return document.querySelector('#summary-table tbody tr [data-statistic=total]').textContent"),"37.5");
    assert_eq!(browser.js("return document.querySelector('#summary-table tbody tr [data-statistic=mean]').textContent"),"18.75");
    assert_eq!(browser.js("return document.querySelectorAll('#summary-table tbody tr')[1].querySelector('[data-statistic=missing]').textContent"),"2");
    assert_eq!(browser.js("return document.querySelectorAll('#summary-table tbody tr')[1].querySelector('[data-statistic=total]').textContent"),"—");
    browser.js("document.getElementById('transcript-metric').value='abundance';document.getElementById('transcript-metric').dispatchEvent(new Event('change'));document.getElementById('summary-threshold').value='1.25'");
    browser.click("#compute-summary");
    browser.wait("return !document.getElementById('download-summary').disabled");
    assert_eq!(browser.js("return document.querySelector('#summary-table tbody tr [data-statistic=total]').textContent"),"3.75");
    assert_eq!(browser.js("return document.querySelector('#summary-table tbody tr [data-statistic=above]').textContent"),"1","threshold comparison must be strict");
    browser.js("document.getElementById('row-sort').value='reference';document.getElementById('column-sort').value='reference';document.getElementById('hide-empty').checked=false");
    browser.click("#select-datasets");
    browser.click("#draw");
    browser.wait("return !document.getElementById('download').disabled");
    assert_eq!(
        browser.js("return document.querySelector('#heatmap td').textContent"),
        "1.25"
    );
    assert_eq!(browser.js("return document.querySelectorAll('#heatmap tbody tr')[1].querySelector('td').textContent"),"—");
    let before = browser.js("return window.queries.length");
    browser.js("document.getElementById('value-mode').value='fold';document.getElementById('value-mode').dispatchEvent(new Event('change'))");
    assert_eq!(
        browser.js("return window.queries.length"),
        before,
        "choosing FC without cached means must not start a scan"
    );
    assert_eq!(
        browser.js("return document.getElementById('download').disabled"),
        true
    );
    browser.click("#draw");
    browser.wait("return !document.getElementById('download').disabled");
    assert_eq!(
        browser.js("return document.querySelector('#heatmap td').textContent"),
        "0.666667"
    );
    assert_eq!(browser.js("return document.querySelectorAll('#heatmap tbody tr')[1].querySelector('td').textContent"),"—");
    assert_eq!(browser.js("return document.querySelectorAll('#heatmap tbody tr')[2].querySelector('td').textContent"),"0");
    browser.click("#compute-summary");
    browser.wait("return !document.getElementById('download-summary').disabled && document.getElementById('cancel').hidden");
    assert_eq!(browser.js("return document.querySelector('#summary-table tbody tr [data-statistic=total]').textContent"),"3.75","summary metrics are not fold transformed");
    browser.js("document.getElementById('value-mode').value='log-fold';document.getElementById('value-mode').dispatchEvent(new Event('change'))");
    assert_eq!(browser.js("return document.querySelectorAll('#heatmap tbody tr')[2].querySelector('td').textContent"),"−∞");
    browser.click("#ranking-tab");
    assert_eq!(
        browser.js("return document.getElementById('ranking-values').disabled"),
        true
    );
    browser.click("#draw-ranking");
    browser.wait("return !document.getElementById('download-ranking').disabled");
    assert_eq!(
        browser.js("return document.getElementById('ranking-title').textContent"),
        "Most expressed transcripts · A"
    );
    assert_eq!(browser.js("return document.getElementById('ranking-details').textContent.includes('2 known, 1 missing')"),true);
    assert_eq!(browser.js("return document.querySelector('#ranking-chart g[data-target=\"0\"]').textContent.includes('1.25')"),true);
    browser.js("document.getElementById('ranking-transcript-metric').value='assigned_observations';document.getElementById('ranking-transcript-metric').dispatchEvent(new Event('change'))");
    browser.click("#draw-ranking");
    browser.wait("return !document.getElementById('download-ranking').disabled");
    assert_eq!(browser.js("return document.querySelector('#ranking-chart g[data-target=\"0\"]').textContent.includes('12.5')"),true);
    browser.click("#target-focus-tab");
    browser.wait("return !document.getElementById('download-target-focus').disabled && document.getElementById('cancel').hidden");
    assert_eq!(browser.js("return document.querySelector('#target-focus-cards [data-statistic=total] strong').title"),"37.5");
    browser.js("document.getElementById('target-focus-name').value='T2';document.getElementById('target-focus-name').dispatchEvent(new Event('change'))");
    browser.click("#load-target-focus");
    browser.wait("return !document.getElementById('download-target-focus').disabled && document.getElementById('cancel').hidden");
    assert_eq!(browser.js("return document.querySelector('#target-focus-cards [data-statistic=known] strong').textContent"),"0");
    assert_eq!(browser.js("return document.querySelector('#target-focus-cards [data-statistic=missing] strong').textContent"),"2");
    assert_eq!(browser.js("return document.querySelector('#target-focus-cards [data-statistic=total] strong').textContent"),"—");
    assert_eq!(browser.js("return [...document.querySelectorAll('#target-focus-table td')].every(td=>td.textContent==='—')"),true);
    browser.click("#overview-tab");
    browser.click("#compute-overview");
    browser.wait("return !document.getElementById('download-overview').disabled && document.getElementById('cancel').hidden");
    assert_eq!(
        browser.js(
            "return document.querySelector('#overview-cards [data-statistic=total] strong').title"
        ),
        "52.5"
    );
    assert_eq!(browser.js("return document.querySelector('#overview-cards [data-statistic=known] strong').textContent"),"4");
    assert_eq!(browser.js("return document.querySelector('#overview-cards [data-statistic=missing] strong').textContent"),"2");
}
