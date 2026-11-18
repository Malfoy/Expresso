// Summaries retain per-target histograms/samples, never a dataset-by-target matrix.
const compare=(a,b)=>a<b?-1:a>b?1:0;
export function parsePercentiles(text) {
  if (!text.trim()) return [];
  const values=text.split(/[\s,;]+/).filter(Boolean).map(value=>{
    if (!/^\d+(\.\d{1,6})?$/.test(value)) throw new Error('Percentiles must be numbers between 0 and 100, with at most 6 decimal places.');
    const number=Number(value);
    if (number<0 || number>100) throw new Error('Percentiles must be between 0 and 100.');
    return number;
  });
  if (values.length>32) throw new Error('Choose at most 32 percentiles.');
  return [...new Set(values)].sort((a,b)=>a-b);
}
function ratio(n,d) {
  if (d===0n) return null;
  const scaled=(n*1000000n+d/2n)/d;
  const whole=scaled/1000000n, fraction=(scaled%1000000n).toString().padStart(6,'0').replace(/0+$/,'');
  return whole.toString()+(fraction?'.'+fraction:'');
}
function midpoint(a,b) { return typeof a==='bigint' ? ratio(a+b,2n) : (a+b)/2; }
function interpolate(a,b,weight,denominator) {
  return typeof a==='bigint' ? ratio(a*BigInt(denominator-weight)+b*BigInt(weight),BigInt(denominator)) : a+(b-a)*weight/denominator;
}
export class Distribution {
  constructor({integer=true, threshold=integer?0n:0, excludeZeros=false, quantiles=true, method='exact', sampleSize=4096, seed=1}={}) {
    Object.assign(this,{integer,threshold,excludeZeros,quantiles,method,sampleSize});
    this.total=integer?0n:0; this.known=0; this.missing=0; this.above=0; this.nonzero=0; this.n=0;
    this.min=null;this.max=null;this.baseline=null;this.center=0;this.m2=0;
    this.histogram=new Map();this.sample=[];this.seed=seed>>>0;
  }
  random() {
    this.seed=(this.seed+0x6D2B79F5)>>>0;
    let n=this.seed;n=Math.imul(n^(n>>>15),n|1);n^=n+Math.imul(n^(n>>>7),n|61);
    return ((n^(n>>>14))>>>0)/4294967296;
  }
  add(value) {
    if (value===null) { this.missing++;return; }
    if ((this.integer && typeof value!=='bigint') || (!this.integer && (typeof value!=='number' || !Number.isFinite(value))) || value<0)
      throw new Error('Invalid summary abundance.');
    this.known++;this.total+=value;
    if (!this.integer && !Number.isFinite(this.total)) throw new Error('Summary total exceeds floating-point range.');
    if (value>this.threshold) this.above++;
    if (value>0) this.nonzero++;
    if (this.excludeZeros && value==0) return;
    this.n++;
    if (this.min===null || value<this.min) this.min=value;
    if (this.max===null || value>this.max) this.max=value;
    if (this.baseline===null) this.baseline=value;
    // Subtract the BigInt baseline first, so small differences near u64::MAX survive.
    const shifted=Number(value-this.baseline),delta=shifted-this.center;
    this.center+=delta/this.n;this.m2+=delta*(shifted-this.center);
    if (!this.quantiles) return;
    if (this.method==='exact') this.histogram.set(value,(this.histogram.get(value)||0)+1);
    else if (this.sample.length<this.sampleSize) this.sample.push(value);
    else { const position=Math.floor(this.random()*this.n);if (position<this.sampleSize) this.sample[position]=value; }
  }
  result(percentiles=[], interpolation='linear') {
    const result={known:this.known,missing:this.missing,nonzero:this.nonzero,above:this.above,
      above_percent:this.known?100*this.above/this.known:null,total:this.known?this.total:null,
      distribution_n:this.n,mean:this.n?(this.integer?ratio(this.total,BigInt(this.n)):this.total/this.n):null,
      min:this.min,max:this.max,stddev:this.n?Math.sqrt(Math.max(0,this.m2/this.n)):null,
      median:null,percentiles:{},quantiles_approximate:this.quantiles && this.method==='sample' && this.n>this.sampleSize,
      sample_n:this.method==='sample'?this.sample.length:this.n};
    if (!this.n || !this.quantiles) return result;
    const entries=this.method==='exact'?[...this.histogram.entries()].sort((a,b)=>compare(a[0],b[0])):this.sample.slice().sort(compare).map(v=>[v,1]);
    let cumulative=0;
    const ends=entries.map(([_,count])=>cumulative+=count),n=cumulative;
    const at=position=>{
      let low=0,high=ends.length-1;
      while (low<high) {const middle=(low+high)>>1;if (ends[middle]>position) high=middle;else low=middle+1;}
      return entries[low][0];
    };
    result.median=n%2?at(Math.floor(n/2)):midpoint(at(n/2-1),at(n/2));
    for (const percentile of percentiles) {
      // P0/P100 always use observed endpoints, including in sampled mode.
      if (percentile===0 || percentile===100) { result.percentiles[percentile]=percentile===0?this.min:this.max;continue; }
      if (interpolation==='nearest-rank') result.percentiles[percentile]=at(Math.max(0,Math.ceil(n*percentile/100)-1));
      else {
        const denominator=100000000,scaled=BigInt(n-1)*BigInt(Math.round(percentile*1000000));
        const low=Number(scaled/BigInt(denominator)),remainder=Number(scaled%BigInt(denominator));
        result.percentiles[percentile]=interpolate(at(low),at(Math.min(n-1,low+1)),remainder,denominator);
      }
    }
    return result;
  }
}

export function summaryChunkSize(datasets,{method='exact',sampleSize=4096,quantiles=true}={}) {
  // At most about one million retained values per batch. With >1M datasets,
  // exact mode necessarily retains the distinct distribution of one target.
  const perTarget=quantiles?(method==='sample'?Math.min(sampleSize,datasets):datasets):1;
  return Math.max(1,Math.min(128,Math.floor(1000000/Math.max(1,perTarget))));
}

export function formatStatistic(value, formatter, mode='raw', format='full') {
  if (value===null || value===undefined) return '—';
  if (typeof value==='string') {
    if (mode==='raw' && format==='full') {
      const [whole,fraction]=value.split('.');return BigInt(whole).toLocaleString()+(fraction?'.'+fraction:'');
    }
    return formatter(Number(value),mode,format);
  }
  return formatter(value,mode,format);
}

export function compareStatistics(a,b) {
  if (a===null || a===undefined) return b===null || b===undefined?0:1;
  if (b===null || b===undefined) return -1;
  if (typeof a==='bigint' && typeof b==='bigint') return compare(a,b);
  if (['string','bigint'].includes(typeof a) && ['string','bigint'].includes(typeof b)) {
    const parts=v=>{const [whole,fraction='']=v.toString().split('.');return [BigInt(whole),fraction.padEnd(6,'0')];};
    const x=parts(a),y=parts(b);return compare(x[0],y[0])||compare(x[1],y[1]);
  }
  return compare(Number(a),Number(b));
}
