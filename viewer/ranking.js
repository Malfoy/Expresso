import {palette} from './display.js';
// A bounded heap keeps only the top N targets while a dataset is scanned in chunks.
const better=(a,b)=>a.value>b.value||(a.value===b.value&&a.id<b.id);
export class Ranking {
  constructor(limit=20,integer=true) {
    if (!Number.isSafeInteger(limit)||limit<1||limit>1000) throw new Error('Top N must be an integer from 1 to 1,000.');
    this.limit=limit;this.integer=integer;this.heap=[];
    this.total=integer?0n:0;this.known=0;this.missing=0;this.expressed=0;this.candidates=0;
  }
  add(id,value,candidate=true) {
    if (value===null) {this.missing++;return;}
    if ((this.integer&&typeof value!=='bigint')||(!this.integer&&(typeof value!=='number'||!Number.isFinite(value)))||value<0)
      throw new Error('Invalid ranking abundance.');
    this.known++;this.total+=value;
    if (!this.integer&&!Number.isFinite(this.total)) throw new Error('Ranking total exceeds floating-point range.');
    if (value<=0) return;
    this.expressed++;
    if (!candidate) return;
    this.candidates++;
    const entry={id,value},heap=this.heap;
    if (heap.length<this.limit) {
      heap.push(entry);let i=heap.length-1;
      while (i>0) {
        const parent=(i-1)>>1;
        if (!better(heap[parent],heap[i])) break;
        [heap[parent],heap[i]]=[heap[i],heap[parent]];i=parent;
      }
    } else if (better(entry,heap[0])) {
      heap[0]=entry;let i=0;
      while (i*2+1<heap.length) {
        let child=i*2+1;
        if (child+1<heap.length&&better(heap[child],heap[child+1])) child++;
        if (!better(heap[i],heap[child])) break;
        [heap[i],heap[child]]=[heap[child],heap[i]];i=child;
      }
    }
  }
  result() {
    const rows=this.heap.slice().sort((a,b)=>better(a,b)?-1:better(b,a)?1:0);
    const displayed=rows.reduce((sum,row)=>sum+row.value,this.integer?0n:0);
    const other=this.integer?this.total-displayed:Math.max(0,this.total-displayed);
    return {rows,total:this.total,displayed,other,known:this.known,missing:this.missing,expressed:this.expressed,candidates:this.candidates};
  }
}
export function share(value,total) {return total>0?Number(value)/Number(total)*100:null;}
export function percentage(value) {
  if (value===null) return '—';
  if (value>0&&value<0.01) return value.toExponential(2)+'%';
  return value.toLocaleString(undefined,{minimumFractionDigits:2,maximumFractionDigits:2})+'%';
}

const namespace='http://www.w3.org/2000/svg';
function element(tag,attrs={},text) {
  const node=document.createElementNS(namespace,tag);
  for (const [key,value] of Object.entries(attrs)) node.setAttribute(key,value);
  if (text!==undefined) node.textContent=text;
  return node;
}
export function barplot({rows,title,subtitle,axis='fit',format,valueLabel,denominatorLabel,accent='#2563eb'}) {
  const colors=palette(accent);
  const width=1100,left=235,right=275,top=100,rowHeight=37,plotWidth=width-left-right;
  const height=top+rows.length*rowHeight+75,svg=element('svg',{xmlns:namespace,viewBox:`0 0 ${width} ${height}`,width,height,role:'img','aria-label':title,'font-family':'system-ui, sans-serif'});
  svg.append(element('title',{},title),element('desc',{},subtitle),element('rect',{width,height,fill:'#ffffff'}));
  svg.append(element('text',{x:24,y:32,fill:'#1b3835','font-size':21,'font-weight':650},title));
  svg.append(element('text',{x:24,y:56,fill:'#628278','font-size':12},subtitle));
  const maximum=Math.max(0,...rows.map(row=>row.percent||0));
  const upper=axis==='fixed'?Math.max(100,Math.ceil(maximum/10)*10):maximum>0?maximum*1.06:100;
  for (let i=0;i<=4;i++) {
    const x=left+plotWidth*i/4;
    svg.append(element('line',{x1:x,x2:x,y1:top-15,y2:top+rows.length*rowHeight,stroke:'#e3ebe3','stroke-width':1}));
    svg.append(element('text',{x,y:top+rows.length*rowHeight+23,'text-anchor':'middle',fill:'#628278','font-size':11},percentage(upper*i/4)));
  }
  rows.forEach((row,i)=>{
    const y=top+i*rowHeight,group=element('g',{'data-target':row.id??'other'});
    const label=row.name.length>30?row.name.slice(0,27)+'…':row.name;
    const text=element('text',{x:left-12,y:y+20,'text-anchor':'end',fill:'#254736','font-size':12},label);
    text.append(element('title',{},row.name));group.append(text);
    const barWidth=row.percent===null?0:Math.max(0,plotWidth*row.percent/upper);
    const bar=element('rect',{x:left,y:y+3,width:barWidth,height:25,rx:4,fill:row.other?colors.css(colors.stops[1]):i===0?colors.main:colors.light,'data-share':row.percent??''});
    const details=`${row.name}: ${format(row.value)} ${valueLabel}; ${percentage(row.percent)} of ${denominatorLabel}`;
    bar.append(element('title',{},details));group.append(bar);
    group.append(element('text',{x:width-right+14,y:y+20,fill:'#254736','font-size':12},`${format(row.value)} · ${percentage(row.percent)}`));
    group.setAttribute('aria-label',details);svg.append(group);
  });
  if (!rows.length) svg.append(element('text',{x:left+plotWidth/2,y:top+12,'text-anchor':'middle',fill:'#628278','font-size':14},'No positive expression among these candidates'));
  svg.append(element('text',{x:left+plotWidth/2,y:height-14,'text-anchor':'middle',fill:'#628278','font-size':12},`Share of ${denominatorLabel} (%)`));
  return svg;
}
