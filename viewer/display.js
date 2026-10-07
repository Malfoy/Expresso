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
  if (mode==='log') {
    const value=Math.log10(Number(count)+1);
    if (format==='scientific') return value.toExponential(2);
    // A log value has fractional information even when full integer counts are selected.
    return value===0 ? '0' : value.toPrecision(3);
  }
  if (format==='human') return human(count);
  if (format==='scientific') {
    const {digits,exponent}=significant(count);
    return `${digits[0]}.${digits.slice(1)}e+${exponent}`;
  }
  return count.toLocaleString();
}

export function needsSummary(mode) {
  return !['reference','name-asc','name-desc'].includes(mode);
}

export function makeSummary(ids) {
  return {index:new Map(ids.map((id,i)=>[id,i])), total:Array(ids.length).fill(0n), max:Array(ids.length).fill(0n), detected:new Uint32Array(ids.length)};
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
