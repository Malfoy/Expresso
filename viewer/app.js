const $ = id => document.getElementById(id);
let file, catalog, targetSelection = new Set(), datasetSelection = new Set(), matrix, worker;
let nextRequest = 0, busy = false, selectionVersion = 0;
const pending = new Map();
function status(message, error = false) { $('status').textContent = message; $('status').classList.toggle('error', error); }
function resetWorker() {
  if (worker) worker.terminate();
  for (const {reject} of pending.values()) reject(new Error('Index changed'));
  pending.clear();
  worker = new Worker('./worker.js', {type:'module'});
  worker.onmessage = ({data}) => {
    const task = pending.get(data.id);
    if (!task) return;
    pending.delete(data.id);
    if (data.error) task.reject(new Error(data.error)); else task.resolve(data.values);
  };
  worker.onerror = error => {
    for (const {reject} of pending.values()) reject(new Error(error.message || 'WASM worker failed'));
    pending.clear();
  };
}
function decode(dataset, selected) {
  return new Promise((resolve, reject) => {
    const id = ++nextRequest;
    pending.set(id, {resolve,reject});
    worker.postMessage({id,file,offset:dataset.offset,length:dataset.length,
      request:{targets:catalog.targets.length, reference:catalog.reference, selected}});
  });
}
function invalidate() {
  selectionVersion++;
  matrix = undefined;
  $('download').disabled = true;
  $('heatmap').hidden = true;
  $('legend').hidden = true;
  $('empty').hidden = false;
}
function choices(kind) {
  const targets = kind === 'targets';
  const entries = targets ? catalog.targets : catalog.datasets;
  const selection = targets ? targetSelection : datasetSelection;
  const query = $(kind + '-search').value.trim().toLowerCase();
  const root = $(kind); root.replaceChildren();
  let matched = 0;
  entries.forEach((entry, i) => {
    if (query && !entry.name.toLowerCase().includes(query) && !(targets && '#' + (i+1) === query)) return;
    matched++;
    if (matched > 150) return;
    const label = document.createElement('label');
    const check = document.createElement('input'); check.type = 'checkbox'; check.checked = selection.has(i); check.dataset.index = i;
    check.addEventListener('change', () => {
      if (check.checked) selection.add(i); else selection.delete(i);
      invalidate(); updateNotes();
    });
    const text = document.createElement('span');
    text.textContent = targets ? `#${i+1} ${entry.name}` : entry.name + (entry.global ? ' · global' : entry.kind==='statistic' ? ' · statistic' : '');
    label.append(check,text); root.append(label);
  });
  $(kind + '-note').dataset.matches = matched;
  updateNotes();
}
function updateNotes() {
  for (const [kind, selection] of [['targets',targetSelection],['datasets',datasetSelection]]) {
    const matches = Number($(kind + '-note').dataset.matches || 0);
    $(kind + '-note').textContent = `${selection.size} selected · ${matches.toLocaleString()} matches${matches>150 ? ' (first 150 shown)' : ''}`;
  }
}
function applyList(kind) {
  if (!catalog) return;
  const targets = kind === 'targets';
  const names = new Set($(targets ? 'target-list' : 'dataset-list').value.split(/\r?\n/).map(x=>x.trim()).filter(x=>x && !x.startsWith('//')));
  const selection = targets ? targetSelection : datasetSelection;
  const entries = targets ? catalog.targets : catalog.datasets;
  const selected = new Set(), found = new Set();
  entries.forEach((entry,i) => {
    for (const key of targets ? [entry.name,'#'+(i+1)] : [entry.name]) {
      if (names.has(key)) { selected.add(i); found.add(key); }
    }
  });
  const missing = [...names].filter(x=>!found.has(x));
  if (missing.length) { status('Unknown selections: '+missing.slice(0,5).join(', '),true); return; }
  for (const i of selected) selection.add(i);
  invalidate(); choices(kind); status(`Added ${selected.size} ${kind} to the selection.`);
}
function u64(value, fallback) {
  if (!value.trim()) return fallback;
  if (!/^\d+$/.test(value.trim())) throw new Error('Count bounds must be nonnegative integers.');
  const n = BigInt(value.trim());
  if (n > 18446744073709551615n) throw new Error('Count bound exceeds u64.');
  return n;
}
function bounds() {
  const min = u64($('min').value,0n), max = u64($('max').value,18446744073709551615n);
  if (min>max) throw new Error('Minimum count exceeds maximum count.');
  return {min,max};
}
function color(fraction) {
  const stops = [[239,247,217],[167,206,130],[69,149,107],[23,83,77]];
  const position = Math.max(0,Math.min(1,fraction))*3;
  const start = Math.min(2,Math.floor(position)), ratio = position-start;
  const rgb = stops[start].map((v,i)=>Math.round(v+(stops[start+1][i]-v)*ratio));
  return `rgb(${rgb.join(',')})`;
}
function render() {
  if (!matrix) return;
  const {min,max} = bounds();
  let highest = 0n;
  for (const column of matrix.columns) for (const value of column) if (value>=min && value<=max && value>highest) highest=value;
  const logarithmic = $('scale').value === 'log';
  const transform = value => logarithmic ? Math.log10(Number(value)+1) : Number(value);
  const denominator = transform(highest) || 1;
  const table = document.createElement('table'), head = document.createElement('thead'), header = document.createElement('tr');
  const corner = document.createElement('th'); corner.textContent = catalog.level === 'gene' ? 'Gene' : 'Exon'; header.append(corner);
  for (const i of matrix.datasets) { const th=document.createElement('th'); th.textContent=catalog.datasets[i].name; header.append(th); }
  head.append(header); table.append(head);
  const body=document.createElement('tbody');
  matrix.targets.forEach((id,row) => {
    const tr=document.createElement('tr'), th=document.createElement('th');
    th.scope='row'; th.textContent=catalog.targets[id].name; th.title=`#${id+1}`; tr.append(th);
    matrix.columns.forEach((column,col) => {
      const value=column[row], visible=value>=min && value<=max;
      const td=document.createElement('td'); td.tabIndex=0;
      const fraction=transform(value)/denominator;
      td.style.background=visible ? color(fraction) : '#edf0ed'; td.style.color=visible && fraction>.62 ? '#fff' : '#254736';
      td.textContent=visible ? value.toLocaleString() : '—';
      td.title=`${catalog.targets[id].name} / ${catalog.datasets[matrix.datasets[col]].name}: ${value}${visible ? '' : ' (outside bounds)'}`;
      td.setAttribute('aria-label',td.title); tr.append(td);
    }); body.append(tr);
  });
  table.append(body); $('heatmap').replaceChildren(table);
  $('heatmap').hidden=false; $('empty').hidden=true; $('legend').hidden=false; $('download').disabled=false;
  $('legend-max').textContent=highest.toLocaleString();
  $('view-title').textContent=`${matrix.targets.length} ${catalog.level === 'gene' ? 'genes' : 'exons'} × ${matrix.datasets.length} datasets`;
}
function csvCell(value) { return '"'+String(value).replaceAll('"','""')+'"'; }
$('download').addEventListener('click', () => {
  try {
    if (!matrix) return;
    const {min,max}=bounds(), gene=catalog.level==='gene';
    const rows=[[gene?'gene_id':'exon_id',gene?'gene_name':'exon_name',...matrix.datasets.map(i=>catalog.datasets[i].name)].map(csvCell).join(',')];
    matrix.targets.forEach((i,row) => rows.push([i+1,catalog.targets[i].name,...matrix.columns.map(col=>col[row]>=min && col[row]<=max ? col[row].toString() : '')].map(csvCell).join(',')));
    const url=URL.createObjectURL(new Blob([rows.join('\r\n')+'\r\n'],{type:'text/csv;charset=utf-8'}));
    const a=document.createElement('a'); a.href=url; a.download='expresso-selection.csv'; a.click(); setTimeout(()=>URL.revokeObjectURL(url),1000);
  } catch (error) { status(error.message,true); }
});
$('draw').addEventListener('click', async () => {
  if (!catalog || busy) return;
  try {
    bounds();
    const targets=[...targetSelection].sort((a,b)=>a-b), datasets=[...datasetSelection].sort((a,b)=>a-b);
    if (!targets.length || !datasets.length) throw new Error('Select at least one target and one dataset.');
    if (targets.length>200 || datasets.length>100) throw new Error('Select at most 200 targets and 100 datasets per view. Use CLI export for larger selections.');
    busy=true; $('draw').disabled=true; $('file').disabled=true;
    matrix=undefined; $('download').disabled=true;
    const version=selectionVersion, columns=[];
    for (const [j,i] of datasets.entries()) {
      status(`Decoding dataset ${j+1} of ${datasets.length}: ${catalog.datasets[i].name}`);
      columns.push((await decode(catalog.datasets[i],targets)).map(BigInt));
    }
    // Store only selected values; full dataset vectors are freed in the worker.
    if (version!==selectionVersion) throw new Error('Selection changed during decoding. Visualize the new selection.');
    matrix={targets,datasets,columns}; render(); status('Selection ready. Counts decoded locally with WebAssembly.');
  } catch (error) { status(error.message,true); }
  finally { busy=false; $('draw').disabled=false; $('file').disabled=false; }
});
$('file').addEventListener('change', async event => {
  try {
    const selected=event.target.files[0]; if (!selected) return;
    invalidate(); $('controls').hidden=true; catalog=undefined;
    status('Reading index catalogue…');
    const header=new Uint8Array(await selected.slice(0,24).arrayBuffer());
    if (header.length!==24 || new TextDecoder().decode(header.slice(0,8))!=='EXPREAI1') throw new Error('Expected an EXPRESSO .eai file produced by expresso pack.');
    const view=new DataView(header.buffer), offset=Number(view.getBigUint64(8,true)), length=Number(view.getBigUint64(16,true));
    if (!Number.isSafeInteger(offset) || !Number.isSafeInteger(length) || offset<24 || length<2 || length>268435456 || offset+length+32!==selected.size) throw new Error('Invalid or oversized index catalogue (maximum 256 MiB).');
    const raw=await selected.slice(offset,offset+length).arrayBuffer();
    const expected=new Uint8Array(await selected.slice(offset+length).arrayBuffer());
    const actual=new Uint8Array(await crypto.subtle.digest('SHA-256',raw));
    if (!actual.every((v,i)=>v===expected[i])) throw new Error('Index catalogue checksum mismatch.');
    const data=JSON.parse(new TextDecoder().decode(raw));
    if (data.version!==1 || !['exon','gene'].includes(data.level) || !Array.isArray(data.targets) || !Array.isArray(data.datasets) || !/^[a-f0-9]{64}$/.test(data.reference)) throw new Error('Invalid index catalogue.');
    if (!data.targets.length || !data.datasets.length || data.targets.some(t=>typeof t.name!=='string')) throw new Error('Index has no valid targets/datasets.');
    let end=24;
    for (const entry of data.datasets) {
      if (typeof entry.name!=='string' || typeof entry.global!=='boolean' || !['dataset','global','statistic'].includes(entry.kind) || !Number.isSafeInteger(entry.offset) || !Number.isSafeInteger(entry.length) || entry.length<=0 || entry.offset!==end || entry.offset+entry.length>offset) throw new Error('Invalid dataset offset/length.');
      end=entry.offset+entry.length;
    }
    if (end!==offset) throw new Error('Inconsistent index length.');
    file=selected; catalog=data; targetSelection.clear(); datasetSelection.clear();
    $('targets-search').value=''; $('datasets-search').value=''; resetWorker();
    $('controls').hidden=false; $('target-count').textContent=data.targets.length.toLocaleString();
    $('target-kind').textContent=data.level==='gene'?'GENES':'EXONS'; $('targets-label').textContent=data.level==='gene'?'Genes':'Exons';
    $('dataset-count').textContent=data.datasets.filter(d=>d.kind==='dataset').length.toLocaleString();
    $('file-size').textContent=(selected.size/1e9).toFixed(2)+' GB'; choices('targets'); choices('datasets');
    status(`Loaded ${selected.name}.${data.coverage && data.coverage.complete===false ? ' Incomplete source coverage: failed datasets were excluded.' : ''}`);
  } catch (error) { status(error.message,true); }
});
for (const kind of ['targets','datasets']) $(kind+'-search').addEventListener('input',()=>{if(catalog) choices(kind);});
$('apply-targets').addEventListener('click',()=>applyList('targets'));
$('apply-datasets').addEventListener('click',()=>applyList('datasets'));
for (const id of ['min','max','scale']) $(id).addEventListener('change',()=>{try{render();}catch(error){$('download').disabled=true;status(error.message,true);}});
$('clear').addEventListener('click',()=>{targetSelection.clear();datasetSelection.clear();invalidate();choices('targets');choices('datasets');});
