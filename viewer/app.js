import {formatAbundance, makeSummary, needsSummary, sortIds} from './display.js';
const $ = id => document.getElementById(id);
let file, catalog, targetSelection = new Set(), datasetSelection = new Set(), matrix, worker;
let nextRequest = 0, busy = false, selectionVersion = 0;
const PAGE_ROWS = 100, PAGE_COLUMNS = 25;
let pageRows = 0, pageColumns = 0;
let summaryCache;
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
function invalidate(selectionChanged=true) {
  if (selectionChanged) { selectionVersion++; summaryCache=undefined; }
  matrix = undefined;
  $('download').disabled = true;
  $('heatmap').hidden = true;
  $('legend').hidden = true;
  $('empty').hidden = false;
  $('paging').hidden = true;
  $('export-selection').disabled = true;
  $('view-details').textContent='';
}
function setBusy(value) {
  busy=value;
  for (const id of ['draw','file','hide-empty','row-sort','column-sort','clear','apply-targets','apply-datasets','select-targets','select-datasets','apply-view']) $(id).disabled=value;
  for (const input of document.querySelectorAll('.choices input')) input.disabled=value;
  $('cancel').hidden=!value;
  if (!value) {
    $('download').disabled=!matrix?.columns || $('heatmap').hidden;
    $('export-selection').disabled=!matrix?.columns || $('heatmap').hidden;
  }
}
function checkVersion(version) {
  if (version!==selectionVersion) throw new Error('Operation canceled.');
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
    const check = document.createElement('input'); check.type = 'checkbox'; check.checked = selection.has(i); check.dataset.index = i; check.disabled=busy;
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
  if (!matrix?.columns) return;
  const {min,max} = bounds();
  let highest = 0n;
  for (const column of matrix.columns) for (const value of column) if (value>=min && value<=max && value>highest) highest=value;
  const logarithmic = $('scale').value === 'log';
  const valueMode=$('value-mode').value, numberFormat=$('number-format').value;
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
      td.textContent=visible ? formatAbundance(value,valueMode,numberFormat) : '—';
      td.title=`${catalog.targets[id].name} / ${catalog.datasets[matrix.datasets[col]].name}: ${value}${valueMode==='log' ? ' · log₁₀(count + 1): '+formatAbundance(value,'log','scientific') : ''}${visible ? '' : ' (outside bounds)'}`;
      td.setAttribute('aria-label',td.title); tr.append(td);
    }); body.append(tr);
  });
  table.append(body); $('heatmap').replaceChildren(table);
  $('heatmap').hidden=false; $('empty').hidden=true; $('legend').hidden=false; $('download').disabled=busy; $('export-selection').disabled=busy;
  $('legend-max').textContent=formatAbundance(highest,'raw',numberFormat)+' counts';
  $('view-details').textContent=`${matrix.hiddenRows ? matrix.hiddenRows.toLocaleString()+' all-zero rows hidden across the selected datasets. ' : ''}Cell values: ${valueMode==='log' ? 'log₁₀(count + 1)' : 'abundance counts'}. CSV exports retain decoded counts.`;
  $('view-title').textContent=`${matrix.allTargets.length.toLocaleString()} ${catalog.level === 'gene' ? 'genes' : 'exons'} × ${matrix.allDatasets.length.toLocaleString()} datasets`;
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
async function loadPage() {
  if (!matrix || busy) return;
  const current=matrix, version=selectionVersion;
  try {
    setBusy(true);
    $('download').disabled=true; $('export-selection').disabled=true;
    for (const id of ['prev-rows','next-rows','prev-columns','next-columns']) $(id).disabled=true;
    bounds();
    const targets=current.allTargets.slice(pageRows,pageRows+PAGE_ROWS);
    const datasets=current.allDatasets.slice(pageColumns,pageColumns+PAGE_COLUMNS);
    const columns=[];
    for (const [j,i] of datasets.entries()) {
      checkVersion(version);
      status(`Loading visible dataset ${j+1} of ${datasets.length}: ${catalog.datasets[i].name}`);
      columns.push((await decode(catalog.datasets[i],targets)).map(BigInt));
    }
    checkVersion(version);
    Object.assign(current,{targets,datasets,columns});
    render(); $('paging').hidden=false; $('export-selection').disabled=false;
    $('row-position').value=pageRows+1; $('column-position').value=pageColumns+1;
    $('page-summary').textContent=`Targets ${pageRows+1}–${pageRows+targets.length} of ${current.allTargets.length.toLocaleString()} · datasets ${pageColumns+1}–${pageColumns+datasets.length} of ${current.allDatasets.length.toLocaleString()}`;
    $('prev-rows').disabled=pageRows===0; $('next-rows').disabled=pageRows+PAGE_ROWS>=current.allTargets.length;
    $('prev-columns').disabled=pageColumns===0; $('next-columns').disabled=pageColumns+PAGE_COLUMNS>=current.allDatasets.length;
    $('view-title').textContent=`${current.allTargets.length.toLocaleString()} ${catalog.level==='gene'?'genes':'exons'} × ${current.allDatasets.length.toLocaleString()} datasets`;
    status('Visible counts ready. Other abundance blocks remain on disk until requested.');
  } catch (error) {
    status(error.message,error.message!=='Operation canceled.');
    if (matrix===current) { $('heatmap').hidden=true; $('legend').hidden=true; }
  } finally { setBusy(false); }
}
async function summarize(targets,datasets,version) {
  const rows=makeSummary(targets), columns=makeSummary(datasets);
  const chunkSize=4096;
  for (const [col,dataset] of datasets.entries()) {
    for (let start=0;start<targets.length;start+=chunkSize) {
      checkVersion(version);
      status(`Scanning selection: dataset ${col+1} of ${datasets.length}, targets ${Math.min(start+chunkSize,targets.length).toLocaleString()} of ${targets.length.toLocaleString()}.`);
      const values=await decode(catalog.datasets[dataset],targets.slice(start,start+chunkSize));
      checkVersion(version);
      values.forEach((text,i)=>{
        const value=BigInt(text), row=start+i;
        rows.total[row]+=value; columns.total[col]+=value;
        if (value>rows.max[row]) rows.max[row]=value;
        if (value>columns.max[col]) columns.max[col]=value;
        if (value>0n) { rows.detected[row]++; columns.detected[col]++; }
      });
    }
  }
  return {version,rows,columns};
}
$('draw').addEventListener('click', async () => {
  if (!catalog || busy) return;
  let targets=[...targetSelection].sort((a,b)=>a-b), datasets=[...datasetSelection].sort((a,b)=>a-b);
  if (!targets.length || !datasets.length) { status('Select at least one target and one dataset.',true); return; }
  const version=selectionVersion, rowSort=$('row-sort').value, columnSort=$('column-sort').value;
  const hideEmpty=$('hide-empty').checked, selectedRows=targets.length;
  invalidate(false);
  try {
    setBusy(true); bounds();
    if (hideEmpty || needsSummary(rowSort) || needsSummary(columnSort)) {
      if (summaryCache?.version!==version) summaryCache=await summarize(targets,datasets,version);
      checkVersion(version);
      if (hideEmpty) targets=targets.filter(id=>summaryCache.rows.detected[summaryCache.rows.index.get(id)]>0);
    }
    targets=sortIds(targets,catalog.targets,rowSort,summaryCache?.rows);
    datasets=sortIds(datasets,catalog.datasets,columnSort,summaryCache?.columns);
    if (!targets.length) {
      $('view-details').textContent=`All ${selectedRows.toLocaleString()} selected rows have zero counts across the selected datasets.`;
      status('No nonzero rows. Disable Hide all-zero rows to display them.');
      return;
    }
    pageRows=0; pageColumns=0;
    matrix={allTargets:targets,allDatasets:datasets,hiddenRows:selectedRows-targets.length};
  } catch(error) { status(error.message,error.message!=='Operation canceled.'); }
  finally { setBusy(false); }
  if (matrix) await loadPage();
});
$('cancel').addEventListener('click',()=>{invalidate();status('Canceling operation…');});
$('apply-view').addEventListener('click',()=>$('draw').click());
for (const id of ['hide-empty','row-sort','column-sort']) $(id).addEventListener('change',()=>{
  invalidate(false); status('View options changed. Visualize selection to apply them.');
});
for (const [id,axis,direction] of [['prev-rows','row',-1],['next-rows','row',1],['prev-columns','column',-1],['next-columns','column',1]]) {
  $(id).addEventListener('click', () => {
    if (!matrix || busy) return;
    if (axis==='row') pageRows+=direction*PAGE_ROWS; else pageColumns+=direction*PAGE_COLUMNS;
    loadPage();
  });
}
$('jump-page').addEventListener('click', () => {
  if (!matrix || busy) return;
  const row=Number($('row-position').value), column=Number($('column-position').value);
  if (!Number.isSafeInteger(row) || !Number.isSafeInteger(column) || row<1 || row>matrix.allTargets.length || column<1 || column>matrix.allDatasets.length) {
    status('Enter positions within the selected target and dataset ranges.',true); return;
  }
  pageRows=Math.floor((row-1)/PAGE_ROWS)*PAGE_ROWS;
  pageColumns=Math.floor((column-1)/PAGE_COLUMNS)*PAGE_COLUMNS;
  loadPage();
});
function selectMatches(kind) {
  if (!catalog || busy) return;
  const target=kind==='targets', entries=target?catalog.targets:catalog.datasets;
  const selection=target?targetSelection:datasetSelection, query=$(kind+'-search').value.trim().toLowerCase();
  entries.forEach((entry,i)=>{
    if (!query || entry.name.toLowerCase().includes(query) || (target && '#'+(i+1)===query)) selection.add(i);
  });
  invalidate(); choices(kind);
}
$('select-targets').addEventListener('click',()=>selectMatches('targets'));
$('select-datasets').addEventListener('click',()=>selectMatches('datasets'));

let exportDirectory, exportFiles=[];
$('export-selection').addEventListener('click', async () => {
  if (!matrix || busy) return;
  const current=matrix, version=selectionVersion;
  let stream, name;
  try {
    const {min,max}=bounds();
    setBusy(true); $('export-selection').disabled=true;
    if (!navigator.storage?.getDirectory) throw new Error('This browser does not support streaming CSV to local storage. Use CLI export.');
    const root=await navigator.storage.getDirectory();
    exportDirectory=await root.getDirectoryHandle('expresso-exports',{create:true});
    name=crypto.randomUUID()+'.csv';
    const handle=await exportDirectory.getFileHandle(name,{create:true});
    stream=await handle.createWritable();
    const gene=catalog.level==='gene';
    await stream.write([gene?'gene_id':'exon_id',gene?'gene_name':'exon_name','dataset','abundance'].map(csvCell).join(',')+'\r\n');
    let rows=0;
    for (const [j,dataset] of current.allDatasets.entries()) {
      status(`Exporting dataset ${j+1} of ${current.allDatasets.length}: ${catalog.datasets[dataset].name}`);
      // Dataset-major order lets the worker reuse one packed vector for every target chunk.
      for (let start=0;start<current.allTargets.length;start+=512) {
        checkVersion(version);
        const targets=current.allTargets.slice(start,start+512), values=await decode(catalog.datasets[dataset],targets);
        checkVersion(version);
        let chunk='';
        values.forEach((text,i)=>{
          const value=BigInt(text);
          if (value>=min && value<=max) {
            const target=targets[i];
            chunk+=[target+1,catalog.targets[target].name,catalog.datasets[dataset].name,text].map(csvCell).join(',')+'\r\n';
            rows++;
          }
        });
        if (chunk) await stream.write(chunk);
      }
    }
    checkVersion(version);
    await stream.close(); stream=undefined;
    checkVersion(version);
    exportFiles.push(name); $('clear-export').hidden=false;
    const output=await handle.getFile(), url=URL.createObjectURL(output);
    const a=document.createElement('a'); a.href=url; a.download='expresso-selection.csv'; a.click();
    // Retain the disk-backed file until the user finishes downloading and clears it.
    exportFiles[exportFiles.length-1]={name,url};
    status(`Exported ${rows.toLocaleString()} matching counts. Download started; clear the temporary CSV after it finishes.`);
  } catch (error) {
    if (stream) await stream.abort().catch(()=>{});
    if (name && exportDirectory && !exportFiles.some(file=>(file.name||file)===name)) await exportDirectory.removeEntry(name).catch(()=>{});
    status(error.message,error.message!=='Operation canceled.');
  } finally {
    setBusy(false);
  }
});
$('clear-export').addEventListener('click', async () => {
  try {
    for (const file of exportFiles) {
      if (file.url) URL.revokeObjectURL(file.url);
      await exportDirectory.removeEntry(file.name||file);
    }
    exportFiles=[]; $('clear-export').hidden=true; status('Temporary CSV files removed.');
  } catch(error) {status(error.message,true);}
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
for (const id of ['min','max','scale','value-mode','number-format']) $(id).addEventListener('change',()=>{
  try {render();} catch(error) {$('download').disabled=true;$('export-selection').disabled=true;status(error.message,true);}
});
$('clear').addEventListener('click',()=>{targetSelection.clear();datasetSelection.clear();invalidate();choices('targets');choices('datasets');});
