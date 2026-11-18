// Presentation and ordering keep u64 counts and aggregate totals as BigInt.
const names = new Intl.Collator(undefined, {numeric:true, sensitivity:'base'});
const compare = (a,b) => a<b ? -1 : a>b ? 1 : 0;

function significant(value) {
  if (value === 0n) return {digits:'000', exponent:0};
  let exponent = value.toString().length-1;
  const shift = exponent-2;
  let rounded;
  if (shift>0) {
    const divisor=10n**BigInt(shift);
    rounded=(value+divisor/2n)/divisor;
  } else rounded=value*10n**BigInt(-shift);
  if (rounded>=1000n) { rounded/=10n; exponent++; }
  return {digits:rounded.toString().padStart(3,'0'),exponent};
}

function human(value) {
  if (value<1000n) return value.toLocaleString();
  const {digits,exponent}=significant(value);
  const group=Math.min(4,Math.floor(exponent/3)), places=exponent-group*3+1;
  const mantissa=places>=3 ? (digits+'0'.repeat(places-3)).replace(/\B(?=(\d{3})+(?!\d))/g,',') : digits.slice(0,places)+'.'+digits.slice(places);
  return mantissa+['','K','M','G','T'][group];
}

export function formatAbundance(count, mode='raw', format='full') {
  if (count===null) return '—';
  if (count===-Infinity) return '−∞';
  if (count===Infinity) return '+∞';
  if (mode==='log') {
    const value=Math.log10(Number(count)+1);
    if (format==='scientific') return value.toExponential(2);
    // A log value has fractional information even when full integer counts are selected.
    return value===0 ? '0' : value.toPrecision(3);
  }
  if (typeof count === 'number') {
    if (format==='scientific') return count.toExponential(2);
    if (format==='human' && count>=1000) {
      const group=Math.min(4,Math.floor(Math.log10(count)/3));
      return Number((count/1000**group).toPrecision(3)).toString()+['','K','M','G','T'][group];
    }
    // Preserve fractional information in normalized values.
    return count.toLocaleString(undefined,{maximumSignificantDigits:6});
  }
  if (format==='human') return human(count);
  if (format==='scientific') {
    const {digits,exponent}=significant(count);
    return `${digits[0]}.${digits.slice(1)}e+${exponent}`;
  }
  return count.toLocaleString();
}

export const isFoldMode=mode=>mode==='fold'||mode==='log-fold';

export function palette(hex='#2563eb') {
  if (!/^#[0-9a-f]{6}$/i.test(hex)) throw new Error('Choose a valid color.');
  const rgb=[1,3,5].map(offset=>parseInt(hex.slice(offset,offset+2),16));
  const blend=(other,weight)=>rgb.map((value,i)=>Math.round(value*(1-weight)+other[i]*weight));
  const css=color=>`rgb(${color.join(',')})`;
  return {rgb,stops:[blend([255,255,255],.94),blend([255,255,255],.65),rgb,blend([0,0,0],.28)],
    main:hex,light:css(blend([255,255,255],.35)),soft:css(blend([255,255,255],.93)),
    border:css(blend([255,255,255],.84)),hover:css(blend([0,0,0],.2)),
    opposite:rgb.map(value=>255-value),css};
}

// Known zeros contribute to the reference mean; missing estimates do not.
// BigInt totals stay exact until the dimensionless ratio is calculated.
export function foldChange(value,total,known,mode='fold') {
  if (value===null || !known || total<=0) return null;
  const ratio=Number(value)/Number(total)*known;
  if (mode!=='log-fold') return ratio;
  if (value===0 || value===0n) return -Infinity;
  // Avoid an underflowed ratio for very small fractional transcript estimates.
  return ratio>0?Math.log2(ratio):Math.log2(Number(value))+Math.log2(known)-Math.log2(Number(total));
}

export function needsSummary(mode) {
  return !['reference','name-asc','name-desc'].includes(mode);
}

export function makeSummary(ids, normalized=false) {
  return {index:new Map(ids.map((id,i)=>[id,i])), total:Array(ids.length).fill(normalized ? 0 : 0n), max:Array(ids.length).fill(normalized ? 0 : 0n), detected:new Uint32Array(ids.length)};
}

const totals = new WeakMap();
export function datasetTotal(dataset) {
  if (totals.has(dataset)) return totals.get(dataset);
  const total=dataset.total_kmers;
  if (!total || typeof total.numerator!=='string' || !/^\d+$/.test(total.numerator) || BigInt(total.numerator)>(1n<<128n)-1n || ![1,1000000].includes(total.scale))
    throw new Error('Dataset k-mer totals are unavailable. Requantify with the current EXPRESSO version.');
  const value=Number(total.numerator)/total.scale;
  totals.set(dataset,value);
  return value;
}

export function normalize(count, target, dataset, mode) {
  if (mode==='none') return count;
  if (dataset.kind==='statistic') throw new Error('Normalization applies to datasets and global sum. Remove raw statistic vectors from the selection.');
  let value=Number(count);
  if (mode==='coverage' || mode==='length-cpm') {
    if (!Number.isSafeInteger(target.kmer_length) || target.kmer_length<0)
      throw new Error('Target k-mer lengths are unavailable. Rebuild the reference index and requantify.');
    if (target.kmer_length===0 && count!==0n) throw new Error('Nonzero abundance has zero exonic k-mer length.');
    value=target.kmer_length===0 ? 0 : value/target.kmer_length;
  }
  if (mode==='cpm' || mode==='length-cpm') {
    const total=datasetTotal(dataset);
    if (total===0 && count!==0n) throw new Error('Nonzero abundance has zero dataset k-mer total.');
    value=total===0 ? 0 : value/total*1000000;
  }
  return value;
}

export function sortIds(ids, entries, mode, summary) {
  return ids.slice().sort((a,b)=>{
    let order=0;
    if (mode.startsWith('name-')) order=names.compare(entries[a].name,entries[b].name);
    else if (needsSummary(mode)) {
      const field=mode.split('-')[0];
      order=compare(summary[field][summary.index.get(a)],summary[field][summary.index.get(b)]);
    }
    if (mode.endsWith('-desc')) order=-order;
    return order || a-b;
  });
}
