import {formatAbundance, makeSummary, needsSummary, sortIds, normalize, datasetTotal, isFoldMode, foldChange, palette} from './display.js';
import {Distribution, parsePercentiles, summaryChunkSize, formatStatistic, compareStatistics} from './statistics.js';
import {Ranking, share, percentage, barplot} from './ranking.js';
import {createExplorer} from './explorer.js';
import {createProgress} from './progress.js';
const $ = id => document.getElementById(id);
let file, catalog, targetSelection = new Set(), datasetSelection = new Set(), matrix, worker;
let nextRequest = 0, busy = false, selectionVersion = 0;
const PAGE_ROWS = 100, PAGE_COLUMNS = 25;
let pageRows = 0, pageColumns = 0;
let summaryCache;
let targetStatistics, statisticsPage=0;
let datasetRanking,rankingPool,rankingPosition=0;
let explorer,activeView='overview';
let themeCache;
let heatmapApplied,statisticsApplied,rankingApplied;
const progress=createProgress($('operation-progress'),$('progress-bar'),$('progress-label'));
const controlKey=ids=>JSON.stringify(ids.map(id=>$(id).type==='checkbox'?$(id).checked:$(id).value));
const heatmapRequest=()=>selectionVersion+'|'+controlKey(['normalization','value-mode','hide-empty','row-sort','column-sort']);
const statisticsRequest=()=>selectionVersion+'|'+controlKey(['summary-scope','summary-values','summary-threshold','summary-percentiles','summary-quantiles','summary-sample-size','summary-interpolation','summary-exclude-zeros'])+'|'+JSON.stringify([...document.querySelectorAll('.summary-metrics input')].map(input=>input.checked));
const rankingRequest=()=>selectionVersion+'|'+controlKey(['ranking-datasets','ranking-position','ranking-limit','ranking-targets','ranking-values','ranking-denominator','ranking-threshold','ranking-percentiles','ranking-quantiles','ranking-transcript-metric']);
function actionState(id,needed) {
  const enabled=!!catalog&&!busy&&needed;
  $(id).disabled=!enabled;$(id).classList.toggle('needs-update',enabled);
}
function updateActions() {
  for (const id of ['draw','apply-view']) actionState(id,heatmapApplied!==heatmapRequest());
  actionState('compute-summary',statisticsApplied!==statisticsRequest());
  actionState('draw-ranking',rankingApplied!==rankingRequest());
  explorer?.updateActions();
}
const colors=()=>themeCache?.main===$('main-color').value?themeCache:(themeCache=palette($('main-color').value));
const targetKind=()=>catalog?.level==='transcript'?'Transcript':catalog?.level==='gene'?'Gene':'Exon';
const metric=()=>catalog?.level==='transcript'?$('transcript-metric').value:'count';
const parseValue=text=>text===null?null:catalog.encoding==='efv'?Number(text):BigInt(text);
const metricLabel=()=>metric()==='abundance'?'EM abundance (coverage rate)':metric()==='assigned_observations'?'assigned k-mer observations':'decoded k-mer counts';
const scalar=()=>catalog.encoding==='efv'?0:0n;
const pending = new Map();
function status(message, error = false) {
  if (busy&&!error) progress.message(message);
  $('status').textContent=busy&&!error?'':message;$('status').classList.toggle('error',error);
}
function resetWorker() {
  if (worker) worker.terminate();
  for (const {reject} of pending.values()) reject(new Error('Index changed'));
  pending.clear();
  worker = new Worker('./worker.js', {type:'module'});
  worker.onmessage = ({data}) => {
    const task = pending.get(data.id);
    if (!task) return;
    pending.delete(data.id);
    if (data.error) task.reject(new Error(data.error)); else {progress.advance();task.resolve(data.values);}
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
    const block=dataset.metrics?dataset.metrics[metric()]:dataset;
    worker.postMessage({id,file,offset:block.offset,length:block.length,
      request:{targets:catalog.targets.length, reference:catalog.reference, selected, encoding:catalog.encoding||'eab'}});
  });
}
function invalidate(selectionChanged=true) {
  heatmapApplied=undefined;
  if (selectionChanged) {
    selectionVersion++; summaryCache=undefined;targetStatistics=undefined;statisticsApplied=undefined;rankingApplied=undefined;
    $('summary-table').hidden=true;$('summary-paging').hidden=true;$('download-summary').disabled=true;
    $('summary-details').textContent='';
    rankingPool=undefined;clearRanking('');updateRankingControls();
    explorer?.selectionChanged();
  }
  matrix = undefined;
  $('download').disabled = true;
  $('heatmap').hidden = true;
  $('legend').hidden = true;
  $('empty').hidden = false;
  $('paging').hidden = true;
  $('export-selection').disabled = true;
  $('view-details').textContent='';
  updateActions();
}
function setBusy(value) {
  if (value!==busy) {if (value) progress.start();else progress.stop();}
  busy=value;
  for (const id of ['draw','file','normalization','transcript-metric','value-mode','hide-empty','row-sort','column-sort','clear','apply-targets','apply-datasets','select-targets','select-datasets','apply-view','compute-summary','download-summary','summary-scope','summary-values','summary-threshold','summary-percentiles','summary-quantiles','summary-sample-size','summary-interpolation','summary-exclude-zeros','summary-sort','summary-direction','summary-color','summary-prev','summary-next']) $(id).disabled=value;
  for (const input of document.querySelectorAll('.summary-metrics input')) input.disabled=value;
  for (const input of document.querySelectorAll('.choices input')) input.disabled=value;
  for (const id of ['comparison-tab','ranking-tab','target-summary-tab','ranking-datasets','ranking-limit','ranking-targets','ranking-values','ranking-denominator','ranking-axis','ranking-format','ranking-other','ranking-prev','ranking-next','ranking-position','ranking-name','ranking-find','ranking-transcript-metric','ranking-threshold','ranking-percentiles','ranking-quantiles','draw-ranking','download-dataset-summary','download-ranking','download-ranking-svg']) $(id).disabled=value;
  for (const element of document.querySelectorAll('.view-tabs button,[role=tabpanel] button,[role=tabpanel] input,[role=tabpanel] select')) element.disabled=value;
  $('cancel').hidden=!value;
  if (!value) {
    $('download').disabled=!matrix?.columns || $('heatmap').hidden;
    $('export-selection').disabled=!matrix?.columns || $('heatmap').hidden;
    $('download-summary').disabled=!targetStatistics;
    $('normalization').disabled=catalog?.level==='transcript';
    $('scale').disabled=isFoldMode(matrix?.valueMode);
    $('summary-sample-size').disabled=$('summary-quantiles').value!=='sample';
    updateRankingControls();
    explorer?.render();
  }
  updateActions();
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
  selection.clear();for (const i of selected) selection.add(i);
  invalidate(); choices(kind); status('');
}
function u64(value, fallback) {
  if (!value.trim()) return fallback;
  if (!/^\d+$/.test(value.trim())) throw new Error('Count bounds must be nonnegative integers.');
  const n = BigInt(value.trim());
  if (n > 18446744073709551615n) throw new Error('Count bound exceeds u64.');
  return n;
}
function bounds() {
  if (catalog.encoding==='efv') {
    const number=(text,fallback)=>{
      if (!text.trim()) return fallback;
      const n=Number(text);if (!Number.isFinite(n) || n<0) throw new Error('Abundance bounds must be finite nonnegative numbers.');return n;
    };
    const min=number($('min').value,0),max=number($('max').value,Infinity);
    if (min>max) throw new Error('Minimum abundance exceeds maximum abundance.');return {min,max};
  }
  const min = u64($('min').value,0n), max = u64($('max').value,18446744073709551615n);
  if (min>max) throw new Error('Minimum count exceeds maximum count.');
  return {min,max};
}
function color(fraction) {
  const stops = colors().stops;
  const position = Math.max(0,Math.min(1,fraction))*3;
  const start = Math.min(2,Math.floor(position)), ratio = position-start;
  const rgb = stops[start].map((v,i)=>Math.round(v+(stops[start+1][i]-v)*ratio));
  return `rgb(${rgb.join(',')})`;
}
function foldColor(value,extent) {
  const position=Number.isFinite(value)?Math.max(-1,Math.min(1,value/extent)):value<0?-1:1;
  const theme=colors(),neutral=[247,247,247],end=position<0?theme.opposite:theme.rgb;
  return `rgb(${neutral.map((v,i)=>Math.round(v+(end[i]-v)*Math.abs(position))).join(',')})`;
}
function foreground(background) {
  const rgb=background.startsWith('#')?palette(background).rgb:background.match(/\d+/g).map(Number);
  const linear=rgb.map(value=>{
    const channel=value/255;return channel<=.04045?channel/12.92:((channel+.055)/1.055)**2.4;
  });
  return .2126*linear[0]+.7152*linear[1]+.0722*linear[2]<=.183?'#fff':'#172b4d';
}
function referenceFor(current,id) {
  const summary=current.reference,row=summary?.index.get(id);
  return row===undefined?{total:scalar(),known:0}:{total:summary.referenceTotal[row],known:summary.referenceKnown[row]};
}
function matrixValue(current,raw,target,dataset,valueMode=current.valueMode) {
  if (raw===null) return null;
  const value=normalize(raw,catalog.targets[target],catalog.datasets[dataset],current.normalization);
  if (!isFoldMode(valueMode)) return value;
  if (catalog.datasets[dataset].kind!=='dataset') return null;
  const reference=referenceFor(current,target);
  return foldChange(value,reference.total,reference.known,valueMode);
}
function render() {
  if (!matrix?.columns) return;
  const {min,max} = bounds();
  const mode=matrix.normalization,valueMode=matrix.valueMode,fold=isFoldMode(valueMode);
  const abundance=(value,row,col)=>matrixValue(matrix,value,matrix.targets[row],matrix.datasets[col]);
  let extent=1;
  let highest = mode==='none' ? scalar() : 0;
  matrix.columns.forEach((column,col)=>column.forEach((value,row)=>{
    if (value!==null && value>=min && value<=max) {
      const normalized=abundance(value,row,col);
      if (normalized!==null && normalized>highest) highest=normalized;
      if (fold) {
        const log=matrixValue(matrix,value,matrix.targets[row],matrix.datasets[col],'log-fold');
        if (log!==null && Number.isFinite(log)) extent=Math.max(extent,Math.abs(log));
      }
    }
  }));
  const logarithmic = $('scale').value === 'log';
  const numberFormat=$('number-format').value;
  const transform = value => logarithmic ? Math.log10(Number(value)+1) : Number(value);
  const denominator = transform(highest) || 1;
  const table = document.createElement('table'), head = document.createElement('thead'), header = document.createElement('tr');
  const corner = document.createElement('th'); corner.textContent = targetKind(); header.append(corner);
  for (const i of matrix.datasets) { const th=document.createElement('th'); th.textContent=catalog.datasets[i].name; header.append(th); }
  head.append(header); table.append(head);
  const body=document.createElement('tbody');
  matrix.targets.forEach((id,row) => {
    const tr=document.createElement('tr'), th=document.createElement('th');
    th.scope='row'; th.textContent=catalog.targets[id].name; th.title=`#${id+1}`; tr.append(th);
    matrix.columns.forEach((column,col) => {
      const raw=column[row],value=raw===null?null:abundance(raw,row,col),visible=value!==null && raw>=min && raw<=max;
      const td=document.createElement('td'); td.tabIndex=0;
      const logFold=fold&&visible?matrixValue(matrix,raw,id,matrix.datasets[col],'log-fold'):null;
      const fraction=fold?Math.abs(logFold)/extent:transform(value)/denominator;
      const background=visible ? fold?foldColor(logFold,extent):color(fraction) : '#edf0ed';
      td.style.background=background;td.style.color=foreground(background);
      td.textContent=visible ? formatAbundance(value,fold?'raw':valueMode,numberFormat) : '—';
      td.title=raw===null?`${catalog.targets[id].name}: missing / unresolved ${metricLabel()}`:`${catalog.targets[id].name} / ${catalog.datasets[matrix.datasets[col]].name}: ${value}${mode==='none' ? '' : ' · raw count: '+raw+' · '+mode}${valueMode==='log' ? ' · log₁₀(value + 1): '+formatAbundance(value,'log','scientific') : ''}${visible ? '' : ' (outside bounds)'}`;
      if (fold && raw!==null) {
        const reference=referenceFor(matrix,id),mean=reference.known?Number(reference.total)/reference.known:null;
        td.title=`${catalog.targets[id].name} / ${catalog.datasets[matrix.datasets[col]].name}: ${formatAbundance(value)} · ${valueMode==='fold'?'fold change':'log₂ fold change'} · target mean: ${mean??'missing'} (${reference.known} known selected real datasets) · source value: ${normalize(raw,catalog.targets[id],catalog.datasets[matrix.datasets[col]],mode)} · raw count: ${raw}`+
          (catalog.datasets[matrix.datasets[col]].kind!=='dataset'?' · global/statistic vectors excluded':reference.total<=0?' · fold change undefined for zero mean':raw<min||raw>max?' (outside bounds)':'');
      }
      td.setAttribute('aria-label',td.title); tr.append(td);
    }); body.append(tr);
  });
  table.append(body); $('heatmap').replaceChildren(table);
  $('heatmap').hidden=false; $('empty').hidden=true; $('legend').hidden=false; $('download').disabled=busy; $('export-selection').disabled=busy;
  $('legend-min').textContent=fold?'≤ '+formatAbundance(valueMode==='fold'?2**(-extent):-extent,'raw',numberFormat):'0';
  $('legend-max').textContent=fold?formatAbundance(valueMode==='fold'?2**extent:extent,'raw',numberFormat):formatAbundance(highest,'raw',numberFormat)+' '+(mode==='none' ? metricLabel() : mode);
  const theme=colors();
  $('legend').querySelector('.gradient').style.background=fold?`linear-gradient(90deg,${theme.css(theme.opposite)},#f7f7f7,${theme.main})`:`linear-gradient(90deg,${theme.stops.map(theme.css).join(',')})`;
  $('legend-note').textContent=fold?'Left: below mean · white: mean · right: above mean. Gray: missing, undefined, excluded, or outside bounds.': 'Gray cells are missing or outside the selected range.';
  $('scale').disabled=busy||fold;
  $('view-details').textContent=matrix.hiddenRows?`${matrix.hiddenRows.toLocaleString()} empty rows hidden`:'';
  $('view-title').textContent=`${matrix.allTargets.length.toLocaleString()} ${targetKind().toLowerCase()}s × ${matrix.allDatasets.length.toLocaleString()} datasets`;
}
function csvCell(value) { return '"'+String(value).replaceAll('"','""')+'"'; }
$('download').addEventListener('click', () => {
  try {
    if (!matrix) return;
    const {min,max}=bounds(), kind=targetKind().toLowerCase();
    const rows=[[kind+'_id',kind+'_name',...matrix.datasets.map(i=>catalog.datasets[i].name)].map(csvCell).join(',')];
    matrix.targets.forEach((i,row) => rows.push([i+1,catalog.targets[i].name,...matrix.columns.map((col,j)=>col[row]!==null && col[row]>=min && col[row]<=max ? matrixValue(matrix,col[row],i,matrix.datasets[j])?.toString()??'' : '')].map(csvCell).join(',')));
    const url=URL.createObjectURL(new Blob([rows.join('\r\n')+'\r\n'],{type:'text/csv;charset=utf-8'}));
    const a=document.createElement('a'); a.href=url; a.download='expresso-selection'+(isFoldMode(matrix.valueMode)?'-'+matrix.valueMode:'')+'.csv'; a.click(); setTimeout(()=>URL.revokeObjectURL(url),1000);
  } catch (error) { status(error.message,true); }
});
async function loadPage(managed=false) {
  if (!matrix || (busy&&!managed)) return;
  const current=matrix, version=selectionVersion;
  try {
    if (!managed) setBusy(true);
    $('download').disabled=true; $('export-selection').disabled=true;
    for (const id of ['prev-rows','next-rows','prev-columns','next-columns']) $(id).disabled=true;
    bounds();
    const targets=current.allTargets.slice(pageRows,pageRows+PAGE_ROWS);
    const datasets=current.allDatasets.slice(pageColumns,pageColumns+PAGE_COLUMNS);
    if (!managed) progress.plan(datasets.length,'Heatmap');
    const columns=[];
    for (const [j,i] of datasets.entries()) {
      checkVersion(version);
      status(`Loading visible dataset ${j+1} of ${datasets.length}: ${catalog.datasets[i].name}`);
      columns.push((await decode(catalog.datasets[i],targets)).map(parseValue));
    }
    checkVersion(version);
    Object.assign(current,{targets,datasets,columns});
    heatmapApplied=heatmapRequest();
    render(); $('paging').hidden=false; $('export-selection').disabled=false;
    $('row-position').value=pageRows+1; $('column-position').value=pageColumns+1;
    $('page-summary').textContent=`Targets ${pageRows+1}–${pageRows+targets.length} of ${current.allTargets.length.toLocaleString()} · datasets ${pageColumns+1}–${pageColumns+datasets.length} of ${current.allDatasets.length.toLocaleString()}`;
    $('prev-rows').disabled=pageRows===0; $('next-rows').disabled=pageRows+PAGE_ROWS>=current.allTargets.length;
    $('prev-columns').disabled=pageColumns===0; $('next-columns').disabled=pageColumns+PAGE_COLUMNS>=current.allDatasets.length;
    $('view-title').textContent=`${current.allTargets.length.toLocaleString()} ${targetKind().toLowerCase()}s × ${current.allDatasets.length.toLocaleString()} datasets`;
    status('');
  } catch (error) {
    status(error.message,error.message!=='Operation canceled.');
    if (matrix===current) { $('heatmap').hidden=true; $('legend').hidden=true; }
  } finally { if (!managed) setBusy(false); }
}
async function summarize(targets,datasets,version,mode) {
  const rows=makeSummary(targets,mode!=='none'||catalog.encoding==='efv'), columns=makeSummary(datasets,mode!=='none'||catalog.encoding==='efv');
  rows.referenceTotal=rows.total.map(()=>mode!=='none'||catalog.encoding==='efv'?0:0n);
  rows.referenceKnown=new Uint32Array(targets.length);
  const chunkSize=4096;
  for (const [col,dataset] of datasets.entries()) {
    for (let start=0;start<targets.length;start+=chunkSize) {
      checkVersion(version);
      status(`Scanning selection: dataset ${col+1} of ${datasets.length}, targets ${Math.min(start+chunkSize,targets.length).toLocaleString()} of ${targets.length.toLocaleString()}.`);
      const values=await decode(catalog.datasets[dataset],targets.slice(start,start+chunkSize));
      checkVersion(version);
      values.forEach((text,i)=>{
        const raw=parseValue(text), row=start+i;
        if (raw===null) return;
        const value=normalize(raw,catalog.targets[targets[row]],catalog.datasets[dataset],mode);
        rows.total[row]+=value; columns.total[col]+=value;
        if (catalog.datasets[dataset].kind==='dataset') {
          rows.referenceTotal[row]+=value;rows.referenceKnown[row]++;
        }
        if (value>rows.max[row]) rows.max[row]=value;
        if (value>columns.max[col]) columns.max[col]=value;
        if (raw>0) { rows.detected[row]++; columns.detected[col]++; }
      });
    }
  }
  return {version,rows,columns};
}
async function drawHeatmap() {
  if (!catalog || busy) return;
  if (heatmapApplied===heatmapRequest()) return;
  activateTab('heatmap',false);
  let targets=[...targetSelection].sort((a,b)=>a-b), datasets=[...datasetSelection].sort((a,b)=>a-b);
  if (!targets.length || !datasets.length) { status('Select at least one target and one dataset.',true); return; }
  const version=selectionVersion, rowSort=$('row-sort').value, columnSort=$('column-sort').value;
  const mode=$('normalization').value;
  const valueMode=$('value-mode').value,fold=isFoldMode(valueMode);
  const hideEmpty=$('hide-empty').checked, selectedRows=targets.length;
  invalidate(false);
  try {
    setBusy(true); bounds();
    const scanning=(fold||hideEmpty||needsSummary(rowSort)||needsSummary(columnSort))&&summaryCache?.version!==version;
    progress.plan((scanning?datasets.length*Math.ceil(targets.length/4096):0)+Math.min(PAGE_COLUMNS,datasets.length),'Heatmap');
    if (fold&&!datasets.some(id=>catalog.datasets[id].kind==='dataset')) throw new Error('Fold change needs at least one selected real dataset. Global/statistic vectors are excluded from the reference mean.');
    // Validate catalogue scalars without querying abundance blocks.
    for (const dataset of datasets) normalize(scalar(),catalog.targets[targets[0]],catalog.datasets[dataset],mode);
    if (mode==='coverage' || mode==='length-cpm') for (const target of targets)
      normalize(0n,catalog.targets[target],catalog.datasets[datasets[0]],mode);
    if (fold || hideEmpty || needsSummary(rowSort) || needsSummary(columnSort)) {
      if (summaryCache?.version!==version) summaryCache=await summarize(targets,datasets,version,mode);
      checkVersion(version);
      if (hideEmpty) targets=targets.filter(id=>summaryCache.rows.detected[summaryCache.rows.index.get(id)]>0);
    }
    targets=sortIds(targets,catalog.targets,rowSort,summaryCache?.rows);
    datasets=sortIds(datasets,catalog.datasets,columnSort,summaryCache?.columns);
    if (!targets.length) {
      $('view-details').textContent=`All ${selectedRows.toLocaleString()} selected rows have zero counts across the selected datasets.`;
      status('No nonzero rows. Disable Hide all-zero rows to display them.');
      heatmapApplied=heatmapRequest();
      return;
    }
    pageRows=0; pageColumns=0;
    matrix={allTargets:targets,allDatasets:datasets,hiddenRows:selectedRows-targets.length,normalization:mode,valueMode,reference:summaryCache?.rows};
    await loadPage(true);
  } catch(error) { status(error.message,error.message!=='Operation canceled.'); }
  finally { setBusy(false); }
}
$('draw').addEventListener('click',drawHeatmap);
$('cancel').addEventListener('click',()=>{invalidate();status('Canceling operation…');});
$('apply-view').addEventListener('click',()=>$('draw').click());
for (const id of ['hide-empty','row-sort','column-sort']) $(id).addEventListener('change',()=>{
  if (heatmapApplied!==heatmapRequest()) invalidate(false);
});
$('normalization').addEventListener('change',()=>{
  invalidate();explorer?.resetSource();status('');
});
$('transcript-metric').addEventListener('change',()=>{
  $('ranking-transcript-metric').value=$('transcript-metric').value;
  invalidate();status('');
  explorer?.resetSource();
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
  selection.clear();
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
    progress.plan(current.allDatasets.length*Math.ceil(current.allTargets.length/512),'Export CSV');
    if (!navigator.storage?.getDirectory) throw new Error('This browser does not support streaming CSV to local storage. Use CLI export.');
    const root=await navigator.storage.getDirectory();
    exportDirectory=await root.getDirectoryHandle('expresso-exports',{create:true});
    name=crypto.randomUUID()+'.csv';
    const handle=await exportDirectory.getFileHandle(name,{create:true});
    stream=await handle.createWritable();
    const kind=targetKind().toLowerCase();
    const field=isFoldMode(current.valueMode)?current.valueMode==='fold'?'fold_change':'log2_fold_change':catalog.level==='transcript'?metric():'abundance';
    await stream.write([kind+'_id',kind+'_name','dataset',field].map(csvCell).join(',')+'\r\n');
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
          const value=parseValue(text);
          if (value!==null && value>=min && value<=max) {
            const target=targets[i];
            const abundance=matrixValue(current,value,target,dataset);
            if (abundance===null) return;
            chunk+=[target+1,catalog.targets[target].name,catalog.datasets[dataset].name,abundance.toString()].map(csvCell).join(',')+'\r\n';
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
    const a=document.createElement('a'); a.href=url; a.download='expresso-selection'+(isFoldMode(current.valueMode)?'-'+current.valueMode:'')+'.csv'; a.click();
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
    invalidate();const openingVersion=selectionVersion;setBusy(true);
    $('controls').hidden=true;$('controls').dataset.ready='false';catalog=undefined;
    status('Reading index catalogue…');
    const header=new Uint8Array(await selected.slice(0,24).arrayBuffer());
    const magic=new TextDecoder().decode(header.slice(0,8));
    if (header.length!==24 || !['EXPREAI1','EXPREAI2'].includes(magic)) throw new Error('Expected an EXPRESSO .eai file produced by expresso pack.');
    const view=new DataView(header.buffer), offset=Number(view.getBigUint64(8,true)), length=Number(view.getBigUint64(16,true));
    if (!Number.isSafeInteger(offset) || !Number.isSafeInteger(length) || offset<24 || length<2 || length>268435456 || offset+length+32!==selected.size) throw new Error('Invalid or oversized index catalogue (maximum 256 MiB).');
    const raw=await selected.slice(offset,offset+length).arrayBuffer();
    const expected=new Uint8Array(await selected.slice(offset+length).arrayBuffer());
    const actual=new Uint8Array(await crypto.subtle.digest('SHA-256',raw));
    checkVersion(openingVersion);
    if (!actual.every((v,i)=>v===expected[i])) throw new Error('Index catalogue checksum mismatch.');
    const data=JSON.parse(new TextDecoder().decode(raw));
    const transcript=data.version===2 && magic==='EXPREAI2' && data.level==='transcript' && data.encoding==='efv' && JSON.stringify(data.metrics)===JSON.stringify(['abundance','assigned_observations']);
    if (!(transcript || (data.version===1 && magic==='EXPREAI1' && ['exon','gene'].includes(data.level))) || !Array.isArray(data.targets) || !Array.isArray(data.datasets) || !/^[a-f0-9]{64}$/.test(data.reference)) throw new Error('Invalid index catalogue.');
    if (!data.targets.length || !data.datasets.length || data.targets.some(t=>typeof t.name!=='string')) throw new Error('Index has no valid targets/datasets.');
    let end=24;
    for (const entry of data.datasets) {
      if (typeof entry.name!=='string' || typeof entry.global!=='boolean' || !['dataset','global','statistic'].includes(entry.kind) || entry.global!==(entry.kind==='global')) throw new Error('Invalid dataset entry.');
      const blocks=transcript?data.metrics.map(metric=>entry.metrics?.[metric]):[entry];
      for (const block of blocks) {
        if (!block || !Number.isSafeInteger(block.offset) || !Number.isSafeInteger(block.length) || block.length<=0 || block.offset!==end || block.offset+block.length>offset) throw new Error('Invalid dataset offset/length.');
        end=block.offset+block.length;
      }
    }
    if (end!==offset) throw new Error('Inconsistent index length.');
    file=selected; catalog=data; targetSelection.clear(); datasetSelection.clear();
    data.targets.forEach((_,id)=>targetSelection.add(id));
    data.datasets.forEach((entry,id)=>{if (entry.kind==='dataset') datasetSelection.add(id);});
    $('normalization').value='none';
    $('normalization').disabled=transcript;$('transcript-metric-controls').hidden=!transcript;
    $('transcript-metric').value='assigned_observations';
    $('ranking-transcript-metric').value='assigned_observations';$('ranking-transcript-controls').hidden=!transcript;
    $('ranking-values').value='none';$('ranking-denominator').value='signal';rankingPosition=0;rankingPool=undefined;
    $('ranking-datasets').value='all';$('ranking-targets').value='all';$('ranking-limit').value='20';
    $('ranking-threshold').value='0';$('ranking-percentiles').value='10, 90, 95, 99';$('ranking-quantiles').value='sample';$('ranking-other').checked=true;
    $('value-mode').value='raw';$('hide-empty').checked=false;
    $('row-sort').value='reference';$('column-sort').value='reference';$('min').value='';$('max').value='';
    $('target-list').value='';$('dataset-list').value='';
    $('targets-search').value=''; $('datasets-search').value=''; resetWorker();
    $('controls').hidden=false; $('target-count').textContent=data.targets.length.toLocaleString();
    $('target-kind').textContent=targetKind().toUpperCase()+'S'; $('targets-label').textContent=targetKind()+'s';
    $('dataset-count').textContent=data.datasets.filter(d=>d.kind==='dataset').length.toLocaleString();
    $('file-size').textContent=(selected.size/1e9).toFixed(2)+' GB'; choices('targets'); choices('datasets');
    $('ranking-title').textContent='Most expressed '+targetKind().toLowerCase()+'s';updateRankingControls();
    $('target-summary-toggle').textContent=targetKind()+' summaries · across datasets';
    $('target-summary-tab').textContent=targetKind()+' list';$('target-focus-tab').textContent=targetKind()+' focus';
    $('summary-title').textContent=targetKind()+' list';$('target-focus-name').value='#1';
    explorer.resetSource();activateTab('overview',false);
    setBusy(false);
    const version=selectionVersion;
    const summarized=explorer.small()?await explorer.computeOverview():true;
    if (version===selectionVersion&&summarized) status(selected.name+(data.coverage?.complete===false?' · incomplete dataset coverage':''));
  } catch (error) { status(error.message,error.message!=='Operation canceled.'); }
  finally {setBusy(false);$('controls').dataset.ready='true';}
});
for (const kind of ['targets','datasets']) $(kind+'-search').addEventListener('input',()=>{if(catalog) choices(kind);});
$('apply-targets').addEventListener('click',()=>applyList('targets'));
$('apply-datasets').addEventListener('click',()=>applyList('datasets'));
for (const id of ['min','max','scale','number-format']) $(id).addEventListener('change',()=>{
  try {render();renderStatistics();explorer?.render();updateActions();} catch(error) {$('download').disabled=true;$('export-selection').disabled=true;status(error.message,true);}
});
$('value-mode').addEventListener('change',()=>{
  try {
    const mode=$('value-mode').value;
    if (matrix && isFoldMode(mode) && !matrix.reference) {
      invalidate(false);
    } else if (matrix) {matrix.valueMode=mode;render();heatmapApplied=heatmapRequest();}
    renderStatistics();updateActions();
  } catch(error) {$('download').disabled=true;$('export-selection').disabled=true;status(error.message,true);}
});
function applyTheme() {
  const theme=colors(),style=document.documentElement.style;
  for (const [key,value] of [['accent',theme.main],['accent-hover',theme.hover],['soft',theme.soft],['border',theme.border]]) style.setProperty('--'+key,value);
  style.setProperty('--on-accent',foreground(theme.css(theme.rgb)));
  $('color-preset').value=[...$('color-preset').options].some(option=>option.value===theme.main)?theme.main:'custom';
  render();renderRanking();renderStatistics();explorer?.render();
}
$('main-color').addEventListener('input',applyTheme);
$('color-preset').addEventListener('change',()=>{
  if ($('color-preset').value==='custom') {$('main-color').focus();return;}
  $('main-color').value=$('color-preset').value;applyTheme();
});
applyTheme();
$('clear').addEventListener('click',()=>{targetSelection.clear();datasetSelection.clear();invalidate();choices('targets');choices('datasets');});

const statisticLabels={total:'Total abundance',above:'Datasets > X',above_percent:'Datasets > X (%)',
  mean:'Mean',median:'Median',min:'Minimum',max:'Maximum',stddev:'Standard deviation',nonzero:'Nonzero datasets'};
function statisticsOptions() {
  const metrics=[...document.querySelectorAll('.summary-metrics input:checked')].map(input=>input.value);
  if (!metrics.length) throw new Error('Choose at least one statistic.');
  const normalization=$('summary-values').value==='view'?$('normalization').value:'none';
  const integer=catalog.encoding!=='efv' && normalization==='none';
  const text=$('summary-threshold').value.trim();
  const threshold=integer?u64(text,0n):text===''?0:Number(text);
  if (!integer && (!Number.isFinite(threshold) || threshold<0)) throw new Error('Summary threshold must be a finite nonnegative number.');
  const method=$('summary-quantiles').value,sampleSize=Number($('summary-sample-size').value);
  if (method==='sample' && (!Number.isSafeInteger(sampleSize) || sampleSize<256 || sampleSize>65536)) throw new Error('Sample size must be an integer from 256 to 65536.');
  const percentiles=metrics.includes('percentiles')?parsePercentiles($('summary-percentiles').value):[];
  return {metrics,normalization,integer,threshold,method,sampleSize,
    percentiles,quantiles:metrics.includes('median')||percentiles.length>0,
    excludeZeros:$('summary-exclude-zeros').checked,interpolation:$('summary-interpolation').value,
    scope:$('summary-scope').value,metric:metric()};
}
function summaryColumns(options) {
  return [['known','Known datasets'],['missing','Missing datasets'],['distribution_n','Distribution N'],
    ...options.metrics.flatMap(field=>field==='percentiles'?options.percentiles.map(p=>['p'+p,'P'+p]):[[field,statisticLabels[field]]])];
}
function statisticValue(row,field) {return field.startsWith('p')?row.values.percentiles[field.slice(1)]:row.values[field];}
function summaryRows() {
  if (!targetStatistics) return [];
  const field=$('summary-sort').value,direction=$('summary-direction').value;
  return targetStatistics.rows.slice().sort((a,b)=>{
    if (field==='reference') return a.id-b.id;
    if (field==='name') {
      const order=catalog.targets[a.id].name.localeCompare(catalog.targets[b.id].name,undefined,{numeric:true});
      return (direction==='desc'?-order:order)||a.id-b.id;
    }
    const x=statisticValue(a,field),y=statisticValue(b,field);
    // Missing estimates always sort last, in either direction.
    if (x===null || x===undefined) return y===null || y===undefined?a.id-b.id:1;
    if (y===null || y===undefined) return -1;
    const order=compareStatistics(x,y);return (direction==='desc'?-order:order)||a.id-b.id;
  });
}
function renderTargetCatalog() {
  if (!catalog) return;
  const ids=[...targetSelection].sort((a,b)=>a-b);
  statisticsPage=Math.max(0,Math.min(statisticsPage,Math.floor(Math.max(0,ids.length-1)/PAGE_ROWS)*PAGE_ROWS));
  const table=document.createElement('table'),head=document.createElement('thead'),header=document.createElement('tr');
  for (const label of [targetKind(),'Original ID','Summary']) {const th=document.createElement('th');th.textContent=label;header.append(th);}
  head.append(header);table.append(head);const body=document.createElement('tbody');
  for (const id of ids.slice(statisticsPage,statisticsPage+PAGE_ROWS)) {
    const row=document.createElement('tr'),th=document.createElement('th'),link=document.createElement('button');
    link.className='table-link';link.textContent=catalog.targets[id].name;link.disabled=busy;
    link.addEventListener('click',()=>explorer.focusTarget(id));th.append(link);row.append(th);
    for (const text of ['#'+(id+1),'Not computed']) {const td=document.createElement('td');td.textContent=text;row.append(td);}body.append(row);
  }
  table.append(body);$('summary-table').replaceChildren(table);$('summary-table').hidden=false;
  $('summary-paging').hidden=false;$('download-summary').disabled=true;
  $('summary-page-label').textContent=`${targetKind()}s ${ids.length?statisticsPage+1:0}–${Math.min(statisticsPage+PAGE_ROWS,ids.length)} of ${ids.length.toLocaleString()}`;
  $('summary-prev').disabled=busy||statisticsPage===0;$('summary-next').disabled=busy||statisticsPage+PAGE_ROWS>=ids.length;
}
function renderStatistics() {
  if (!targetStatistics) {renderTargetCatalog();return;}
  const {options,datasets}=targetStatistics,columns=summaryColumns(options);
  const select=$('summary-sort'),current=select.value;
  select.replaceChildren(...[['reference','Reference order'],['name','Target name'],...columns].map(([field,label])=>{
    const option=document.createElement('option');option.value=field;option.textContent=label;return option;
  }));
  select.value=[...select.options].some(o=>o.value===current)?current:'reference';
  const rows=summaryRows();
  const colorize=$('summary-color').checked,transform=value=>$('scale').value==='log'?Math.log10(Number(value)+1):Number(value);
  const maxima=new Map(columns.map(([field])=>[field,rows.reduce((max,row)=>{
    const value=statisticValue(row,field);return value===null||value===undefined?max:Math.max(max,transform(value));
  },0)]));
  statisticsPage=Math.max(0,Math.min(statisticsPage,Math.floor((rows.length-1)/PAGE_ROWS)*PAGE_ROWS));
  const page=rows.slice(statisticsPage,statisticsPage+PAGE_ROWS),table=document.createElement('table');
  const header=document.createElement('tr');
  for (const label of [targetKind(),...columns.map(([_,label])=>label)]) {
    const th=document.createElement('th');th.textContent=label;header.append(th);
  }
  const head=document.createElement('thead');head.append(header);table.append(head);
  const body=document.createElement('tbody');
  for (const row of page) {
    const tr=document.createElement('tr'),th=document.createElement('th');th.scope='row';
    const link=document.createElement('button');link.className='table-link';link.textContent=catalog.targets[row.id].name;
    link.disabled=busy;link.addEventListener('click',()=>explorer.focusTarget(row.id));th.append(link);th.title='#'+(row.id+1);tr.append(th);
    for (const [field] of columns) {
      const value=statisticValue(row,field),td=document.createElement('td');
      const abundance=!['known','missing','distribution_n','above','above_percent','nonzero'].includes(field);
      td.textContent=formatStatistic(value,formatAbundance,abundance&&$('value-mode').value==='log'?'log':'raw',$('number-format').value);
      td.title=value===null||value===undefined?'No known abundance values':value.toString();
      td.classList.toggle('missing',value===null||value===undefined);td.dataset.statistic=field;tr.append(td);
      if (colorize && value!==null && value!==undefined && !['known','missing','distribution_n'].includes(field)) {
        const fraction=transform(value)/(maxima.get(field)||1),background=color(fraction);td.style.background=background;td.style.color=foreground(background);
      }
    }
    body.append(tr);
  }
  table.append(body);$('summary-table').replaceChildren(table);$('summary-table').hidden=false;
  $('summary-paging').hidden=false;$('summary-page-label').textContent=`Targets ${statisticsPage+1}–${statisticsPage+page.length} of ${rows.length.toLocaleString()}`;
  $('summary-prev').disabled=busy||statisticsPage===0;$('summary-next').disabled=busy||statisticsPage+PAGE_ROWS>=rows.length;
  const sampled=rows.filter(row=>row.values.quantiles_approximate).length;
  $('summary-details').textContent=`${rows.length.toLocaleString()} targets · ${datasets.length.toLocaleString()} datasets${sampled?' · approximate quantiles':''}`;
  $('download-summary').disabled=busy;
}
async function computeTargetSummaries() {
  if (!catalog || busy) return;
  if (statisticsApplied===statisticsRequest()) return;
  const targets=[...targetSelection].sort((a,b)=>a-b),version=selectionVersion;
  if (!targets.length) {status('Select at least one target to summarize.',true);return;}
  try {
    const options=statisticsOptions();
    const datasets=catalog.datasets.map((entry,id)=>({entry,id})).filter(({entry,id})=>entry.kind==='dataset'&&(options.scope==='all'||datasetSelection.has(id))).map(({id})=>id);
    if (!datasets.length) throw new Error('The summary scope contains no real datasets. Global/statistic vectors are excluded.');
    // Validate normalization before reading any abundance blocks.
    for (const dataset of datasets) normalize(scalar(),catalog.targets[targets[0]],catalog.datasets[dataset],options.normalization);
    if (options.normalization==='coverage'||options.normalization==='length-cpm') for (const id of targets)
      normalize(0n,catalog.targets[id],catalog.datasets[datasets[0]],options.normalization);
    setBusy(true);targetStatistics=undefined;$('summary-table').hidden=true;$('summary-paging').hidden=true;
    const rows=[],chunkSize=summaryChunkSize(datasets.length,options);
    progress.plan(Math.ceil(targets.length/chunkSize)*datasets.length,'Target summaries');
    for (let start=0;start<targets.length;start+=chunkSize) {
      const selected=targets.slice(start,start+chunkSize),distributions=selected.map(id=>new Distribution({...options,seed:id+1}));
      for (const [position,dataset] of datasets.entries()) {
        checkVersion(version);
        status(`Summarizing targets ${start+1}–${start+selected.length} of ${targets.length.toLocaleString()}: dataset ${position+1} of ${datasets.length.toLocaleString()}.`);
        const values=await decode(catalog.datasets[dataset],selected);checkVersion(version);
        if (values.length!==selected.length) throw new Error('Incomplete summary values.');
        values.forEach((text,i)=>{
          const raw=parseValue(text);distributions[i].add(raw===null?null:normalize(raw,catalog.targets[selected[i]],catalog.datasets[dataset],options.normalization));
        });
      }
      for (const [i,id] of selected.entries()) rows.push({id,values:distributions[i].result(options.percentiles,options.interpolation)});
      // The batch's distributions are released before the next batch is decoded.
    }
    checkVersion(version);targetStatistics={rows,options,datasets};statisticsApplied=statisticsRequest();statisticsPage=0;renderStatistics();
    status('');
  } catch(error) {status(error.message,error.message!=='Operation canceled.');}
  finally {setBusy(false);if (targetStatistics) renderStatistics();}
}
$('compute-summary').addEventListener('click',computeTargetSummaries);
function clearStatistics() {
  if (statisticsApplied===statisticsRequest()) {updateActions();return;}
  statisticsApplied=undefined;
  targetStatistics=undefined;$('summary-table').hidden=true;$('summary-paging').hidden=true;$('download-summary').disabled=true;
  $('summary-details').textContent='';
  $('summary-sample-size').disabled=$('summary-quantiles').value!=='sample';
  updateActions();
}
for (const id of ['summary-scope','summary-values','summary-threshold','summary-percentiles','summary-quantiles','summary-sample-size','summary-interpolation','summary-exclude-zeros']) $(id).addEventListener('change',clearStatistics);
for (const input of document.querySelectorAll('.summary-metrics input')) input.addEventListener('change',clearStatistics);
for (const id of ['summary-sort','summary-direction']) $(id).addEventListener('change',()=>{statisticsPage=0;renderStatistics();});
$('summary-color').addEventListener('change',renderStatistics);
$('summary-prev').addEventListener('click',()=>{statisticsPage-=PAGE_ROWS;renderStatistics();});
$('summary-next').addEventListener('click',()=>{statisticsPage+=PAGE_ROWS;renderStatistics();});
$('download-summary').addEventListener('click',()=>{
  if (!targetStatistics || busy) return;
  const kind=targetKind().toLowerCase(),columns=summaryColumns(targetStatistics.options),options=targetStatistics.options;
  const rows=[[kind+'_id',kind+'_name','metric','normalization','dataset_scope','datasets','threshold','exclude_zeros','quantile_method','percentile_definition','quantiles_approximate','quantile_sample_n',...columns.map(([field])=>field)].map(csvCell).join(',')];
  for (const row of summaryRows()) rows.push([row.id+1,catalog.targets[row.id].name,options.metric,options.normalization,options.scope,targetStatistics.datasets.length,options.threshold,options.excludeZeros,options.method,options.interpolation,row.values.quantiles_approximate,row.values.sample_n,
    ...columns.map(([field])=>statisticValue(row,field)??'')].map(csvCell).join(','));
  const url=URL.createObjectURL(new Blob([rows.join('\r\n')+'\r\n'],{type:'text/csv;charset=utf-8'}));
  const link=document.createElement('a');link.href=url;link.download='expresso-target-summaries.csv';link.click();setTimeout(()=>URL.revokeObjectURL(url),10000);
});

function rankingDatasets() {
  if (!catalog) return [];
  if (!rankingPool) rankingPool=catalog.datasets.map((entry,id)=>({entry,id})).filter(({entry,id})=>entry.kind==='dataset'&&($('ranking-datasets').value==='all'||datasetSelection.has(id))).map(({id})=>id);
  return rankingPool;
}
function clearRanking(message='') {
  rankingApplied=undefined;
  datasetRanking=undefined;$('ranking-chart').hidden=true;$('ranking-chart').replaceChildren();
  $('dataset-summary-cards').hidden=true;$('dataset-summary-cards').replaceChildren();
  $('download-dataset-summary').disabled=true;$('download-ranking').disabled=true;$('download-ranking-svg').disabled=true;$('ranking-details').textContent=message;
  updateActions();
}
function updateRankingControls() {
  const pool=rankingDatasets();
  rankingPosition=Math.max(0,Math.min(rankingPosition,pool.length-1));
  $('ranking-position').value=rankingPosition+1;$('ranking-position').max=Math.max(1,pool.length);
  $('ranking-prev').disabled=busy||rankingPosition===0;$('ranking-next').disabled=busy||rankingPosition+1>=pool.length;
  $('download-ranking').disabled=busy||!datasetRanking;$('download-ranking-svg').disabled=busy||!datasetRanking;
  $('download-dataset-summary').disabled=busy||!datasetRanking;
  $('ranking-values').disabled=busy||catalog?.level==='transcript';
  const inputAllowed=catalog?.level!=='transcript'&&$('ranking-values').value==='none';
  $('ranking-denominator').querySelector('option[value=input]').disabled=!inputAllowed;
  if (!inputAllowed) $('ranking-denominator').value='signal';
}
function rankingRows() {
  const report=datasetRanking;
  if (!report) return [];
  const rows=report.rows.map(row=>({...row,name:catalog.targets[row.id].name,percent:share(row.value,report.denominator)}));
  if ($('ranking-other').checked&&report.other>0) rows.push({id:null,name:'Other '+targetKind().toLowerCase()+'s',value:report.other,percent:share(report.other,report.denominator),other:true});
  return rows;
}
function renderRanking() {
  if (!datasetRanking) return;
  const report=datasetRanking,kind=targetKind().toLowerCase(),dataset=catalog.datasets[report.dataset];
  const denominatorLabel=report.denominatorKind==='input'?'full dataset k-mer total':`${catalog.encoding==='efv'?'known ':''}${kind} signal`;
  const valueLabel=report.normalization==='none'?metricLabel():report.normalization;
  const total=report.known?formatAbundance(report.total,'raw',$('ranking-format').value):'—';
  const title=`Most expressed ${kind}s · ${dataset.name}`;
  const subtitle=`Top ${report.rows.length} · ${valueLabel} · total ${kind} signal: ${total} · percentages of ${denominatorLabel}`;
  const svg=barplot({rows:rankingRows(),title,subtitle,axis:$('ranking-axis').value,
    format:value=>formatAbundance(value,'raw',$('ranking-format').value),valueLabel,denominatorLabel,accent:colors().main});
  $('ranking-chart').replaceChildren(svg);$('ranking-chart').hidden=false;
  $('ranking-title').textContent=title;
  $('dataset-summary-title').textContent='Dataset summary · '+dataset.name;
  const stats=report.statistics;
  const fields=[['total','Total expression',stats.total],['nonzero','Expressed '+kind+'s',stats.nonzero],
    ['mean','Mean expression',stats.mean],['median','Median expression',stats.median],
    ['known','Known '+kind+'s',stats.known],['missing','Missing '+kind+'s',stats.missing],
    ['above',`Expression > ${report.threshold}`,stats.above],['above_percent','% above threshold',stats.above_percent],
    ...report.percentiles.map(p=>['p'+p,'P'+p,stats.percentiles[p]])];
  $('dataset-summary-cards').replaceChildren(...fields.map(([field,label,value])=>{
    const card=document.createElement('div'),name=document.createElement('small'),number=document.createElement('strong');
    card.dataset.statistic=field;name.textContent=label;
    number.textContent=field==='above_percent'?percentage(value):formatStatistic(value,formatAbundance,'raw',$('ranking-format').value);
    number.title=value?.toString()??'No known expression';card.append(name,number);return card;
  }));$('dataset-summary-cards').hidden=false;
  const fraction=share(report.displayed,report.denominator);
  $('ranking-details').textContent=`Total: ${total} · top ${report.rows.length}: ${percentage(fraction)} of ${denominatorLabel} · ${report.known.toLocaleString()} known, ${report.missing.toLocaleString()} missing${stats.quantiles_approximate?' · approximate quantiles':''}`;
  $('download-dataset-summary').disabled=busy;$('download-ranking').disabled=busy;$('download-ranking-svg').disabled=busy;
}
const tabs={'overview':['overview-tab','overview-view'],'datasets':['datasets-tab','datasets-view'],
  'targets':['target-summary-tab','target-summary-panel'],'heatmap':['comparison-tab','comparison-view'],
  'dataset-focus':['ranking-tab','ranking-view'],'target-focus':['target-focus-tab','target-focus-view']};
async function activateTab(view,load=true) {
  if (busy&&load) return;
  activeView=view;
  for (const [name,[tab,panel]] of Object.entries(tabs)) {
    const active=name===view;$(panel).hidden=!active;$(tab).setAttribute('aria-selected',String(active));$(tab).tabIndex=active?0:-1;
  }
  if (!catalog||!load) return;
  if (view==='heatmap') {if (!matrix?.columns) await drawHeatmap();else render();}
  else if (view==='dataset-focus') {if (!datasetRanking) await drawRanking();else renderRanking();}
  else if (view==='targets') {renderStatistics();if (!targetStatistics&&explorer.small()) await computeTargetSummaries();}
  else await explorer.enter(view);
}
for (const [name,[tab]] of Object.entries(tabs)) {
  $(tab).addEventListener('click',()=>activateTab(name));
  $(tab).addEventListener('keydown',event=>{
    const names=Object.keys(tabs),index=names.indexOf(name);
    const next=event.key==='ArrowRight'?(index+1)%names.length:event.key==='ArrowLeft'?(index+names.length-1)%names.length:event.key==='Home'?0:event.key==='End'?names.length-1:-1;
    if (next>=0) {event.preventDefault();$(tabs[names[next]][0]).focus();activateTab(names[next]);}
  });
}
async function drawRanking() {
  if (!catalog||busy) return;
  if (rankingApplied===rankingRequest()) return;
  const version=selectionVersion;
  try {
    const pool=rankingDatasets(),position=Number($('ranking-position').value);
    if (!pool.length) throw new Error('Select at least one real dataset, or choose All real datasets. Global/statistic vectors are excluded.');
    if (!Number.isSafeInteger(position)||position<1||position>pool.length) throw new Error('Enter a dataset position within the chosen dataset scope.');
    const datasetId=pool[position-1],dataset=catalog.datasets[datasetId],normalization=$('ranking-values').value;
    const candidateScope=$('ranking-targets').value,candidates=new Set(targetSelection);
    if (candidateScope==='selected'&&!candidates.size) throw new Error('Select at least one target, or rank All targets.');
    const integer=normalization==='none'&&catalog.encoding!=='efv',ranking=new Ranking(Number($('ranking-limit').value),integer);
    const text=$('ranking-threshold').value.trim();
    if (!(integer?/^\d+$/:/^\d+(\.\d+)?$/).test(text)) throw new Error('Dataset summary threshold must be a nonnegative '+(integer?'integer.':'number.'));
    const threshold=integer?BigInt(text):Number(text);
    if (!integer&&!Number.isFinite(threshold)) throw new Error('Dataset summary threshold must be finite.');
    const percentiles=parsePercentiles($('ranking-percentiles').value);
    const distribution=new Distribution({integer,threshold,method:$('ranking-quantiles').value,sampleSize:4096});
    const denominatorKind=$('ranking-denominator').value;
    let denominator;
    if (denominatorKind==='input') {
      if (!integer) throw new Error('Full dataset denominator is available for raw exon/gene counts only.');
      denominator=datasetTotal(dataset);
    }
    normalize(scalar(),catalog.targets[0],dataset,normalization);
    if (normalization==='coverage'||normalization==='length-cpm') for (const target of catalog.targets) normalize(0n,target,dataset,normalization);
    clearRanking('Scanning the chosen dataset…');rankingPosition=position-1;setBusy(true);
    const chunk=4096;
    progress.plan(Math.ceil(catalog.targets.length/chunk),'Dataset summary');
    for (let start=0;start<catalog.targets.length;start+=chunk) {
      checkVersion(version);
      const selected=Array.from({length:Math.min(chunk,catalog.targets.length-start)},(_,i)=>start+i);
      status(`Ranking ${dataset.name}: ${Math.min(start+chunk,catalog.targets.length).toLocaleString()} of ${catalog.targets.length.toLocaleString()} targets.`);
      const values=await decode(dataset,selected);checkVersion(version);
      if (values.length!==selected.length) throw new Error('Incomplete ranking values.');
      values.forEach((text,i)=>{
        const raw=parseValue(text),id=selected[i];
        const value=raw===null?null:normalize(raw,catalog.targets[id],dataset,normalization);
        ranking.add(id,value,candidateScope==='all'||candidates.has(id));distribution.add(value);
      });
    }
    checkVersion(version);
    const result=ranking.result();
    if (denominatorKind==='input'&&denominator===0&&result.total>0) throw new Error('Nonzero gene/target signal has a zero full dataset k-mer total.');
    datasetRanking={...result,dataset:datasetId,normalization,candidateScope,denominatorKind,
      denominator:denominatorKind==='input'?denominator:result.total,metric:metric(),
      statistics:distribution.result(percentiles),threshold,percentiles};
    rankingApplied=rankingRequest();renderRanking();status('');
  } catch(error) {status(error.message,error.message!=='Operation canceled.');}
  finally {setBusy(false);}
}
$('draw-ranking').addEventListener('click',drawRanking);
for (const id of ['ranking-datasets','ranking-limit','ranking-targets','ranking-values','ranking-denominator','ranking-threshold','ranking-percentiles','ranking-quantiles']) $(id).addEventListener('change',()=>{
  if (rankingApplied===rankingRequest()) {updateActions();return;}
  if (id==='ranking-datasets') {rankingPool=undefined;rankingPosition=0;}
  clearRanking();updateRankingControls();
});
for (const id of ['ranking-axis','ranking-format','ranking-other']) $(id).addEventListener('change',renderRanking);
$('ranking-position').addEventListener('change',()=>{if (rankingApplied!==rankingRequest()) clearRanking();});
$('ranking-transcript-metric').addEventListener('change',()=>{
  $('transcript-metric').value=$('ranking-transcript-metric').value;$('transcript-metric').dispatchEvent(new Event('change'));
});
for (const [id,direction] of [['ranking-prev',-1],['ranking-next',1]]) $(id).addEventListener('click',()=>{
  if (busy) return;rankingPosition+=direction;$('ranking-position').value=rankingPosition+1;drawRanking();
});
$('ranking-find').addEventListener('click',()=>{
  if (!catalog||busy) return;
  const name=$('ranking-name').value.trim(),position=rankingDatasets().findIndex(id=>catalog.datasets[id].name===name);
  if (position<0) {status('Dataset not found in the chosen scope: '+name,true);return;}
  rankingPosition=position;$('ranking-position').value=position+1;drawRanking();
});
function downloadRankingBlob(blob,suffix) {
  const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;
  a.download='expresso-'+catalog.datasets[datasetRanking.dataset].name+'-ranking.'+suffix;a.click();setTimeout(()=>URL.revokeObjectURL(url),10000);
}
$('download-dataset-summary').addEventListener('click',()=>{
  if (!datasetRanking||busy) return;
  const report=datasetRanking,stats=report.statistics;
  const fields=['total','known','missing','nonzero','above','above_percent','mean','median','min','max','stddev'];
  const header=['dataset','target_kind','metric','normalization','threshold','quantiles_approximate','quantile_sample_n',...fields,...report.percentiles.map(p=>'p'+p)];
  const row=[catalog.datasets[report.dataset].name,targetKind().toLowerCase(),report.metric,report.normalization,report.threshold,stats.quantiles_approximate,stats.sample_n,
    ...fields.map(field=>stats[field]??''),...report.percentiles.map(p=>stats.percentiles[p]??'')];
  downloadRankingBlob(new Blob([header.map(csvCell).join(',')+'\r\n'+row.map(csvCell).join(',')+'\r\n'],{type:'text/csv;charset=utf-8'}),'statistics.csv');
});
$('download-ranking-svg').addEventListener('click',()=>{
  if (!datasetRanking||busy) return;
  const svg=$('ranking-chart').querySelector('svg');
  downloadRankingBlob(new Blob([new XMLSerializer().serializeToString(svg)],{type:'image/svg+xml;charset=utf-8'}),'svg');
});
$('download-ranking').addEventListener('click',()=>{
  if (!datasetRanking||busy) return;
  const report=datasetRanking,dataset=catalog.datasets[report.dataset],kind=targetKind().toLowerCase();
  const header=['dataset','rank','kind',kind+'_id',kind+'_name','metric','normalization','value','percent_of_total','denominator_kind','denominator','dataset_kmer_numerator','dataset_kmer_scale','all_target_signal','displayed_signal','known_targets','missing_targets'];
  const rows=rankingRows();
  if (!rows.length) rows.push({id:null,name:'Displayed total',value:report.displayed,percent:share(report.displayed,report.denominator),summary:true});
  const csv=[header.map(csvCell).join(',')];
  for (const [i,row] of rows.entries()) csv.push([dataset.name,row.other||row.summary?'':i+1,row.other?'other':row.summary?'summary':'target',row.id===null?'':row.id+1,row.name,
    report.metric,report.normalization,row.value,row.percent??'',report.denominatorKind,report.denominator,
    report.denominatorKind==='input'?dataset.total_kmers.numerator:'',report.denominatorKind==='input'?dataset.total_kmers.scale:'',report.total,report.displayed,report.known,report.missing].map(csvCell).join(','));
  downloadRankingBlob(new Blob([csv.join('\r\n')+'\r\n'],{type:'text/csv;charset=utf-8'}),'csv');
});

explorer=createExplorer({$,controlKey,actionState,progress:(total,label)=>progress.plan(total,label),catalog:()=>catalog,file:()=>file,metric,metricLabel,targetKind,
  normalization:()=>$('normalization').value,selectedTargets:()=>targetSelection,selectedDatasets:()=>datasetSelection,
  parseValue,read:(dataset,selected)=>decode(catalog.datasets[dataset],selected),colors,
  busy:()=>busy,setBusy,status,version:()=>selectionVersion,checkVersion,
  activate:view=>activateTab(view,false),
  focusDataset:async id=>{
    if (busy) return;
    activateTab('dataset-focus',false);$('ranking-datasets').value='all';rankingPool=undefined;
    rankingPosition=rankingDatasets().indexOf(id);$('ranking-values').value=$('normalization').value;
    $('ranking-targets').value='all';clearRanking();updateRankingControls();await drawRanking();
  }
});
activateTab('overview',false);

for (const event of ['input','change']) document.addEventListener(event,updateActions);
updateActions();
