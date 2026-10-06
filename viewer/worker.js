let wasmPromise, cached;
async function module() {
  if (!wasmPromise) wasmPromise = (async () => {
    const response = await fetch('./expresso_viewer.wasm');
    if (!response.ok) throw new Error('WASM decoder missing. Follow the viewer build instructions in the README.');
    return (await WebAssembly.instantiate(await response.arrayBuffer(), {})).instance.exports;
  })();
  return wasmPromise;
}
function call(wasm, request, vector) {
  const bytes = new TextEncoder().encode(JSON.stringify(request));
  let requestPtr, vectorPtr, resultPtr, resultLen;
  try {
    requestPtr = wasm.allocate(bytes.length);
    if (vector) vectorPtr = wasm.allocate(vector.length);
    new Uint8Array(wasm.memory.buffer, requestPtr, bytes.length).set(bytes);
    if (vector) new Uint8Array(wasm.memory.buffer, vectorPtr, vector.length).set(vector);
    const packed = vector ? wasm.prepare_vector(requestPtr, bytes.length, vectorPtr, vector.length) : wasm.select_counts(requestPtr, bytes.length, cached.handle);
    resultPtr = Number(packed & 0xffffffffn);
    resultLen = Number(packed >> 32n);
    const result = JSON.parse(new TextDecoder().decode(new Uint8Array(wasm.memory.buffer, resultPtr, resultLen)));
    if (result.error) throw new Error(result.error);
    return result;
  } finally {
    if (requestPtr !== undefined) wasm.release(requestPtr, bytes.length);
    if (vectorPtr !== undefined) wasm.release(vectorPtr, vector.length);
    if (resultPtr !== undefined) wasm.release(resultPtr, resultLen);
  }
}
async function decode(data) {
  try {
    const wasm = await module();
    const key = `${data.file.size}:${data.file.lastModified}:${data.offset}:${data.length}:${data.request.reference}`;
    if (!cached || cached.key !== key) {
      if (cached) { wasm.release_vector(cached.handle); cached = undefined; }
      const vector = new Uint8Array(await data.file.slice(data.offset, data.offset + data.length).arrayBuffer());
      if (vector.length !== data.length) throw new Error('Truncated dataset block');
      const result = call(wasm, {...data.request, selected:[]}, vector);
      cached = {key, handle:result.handle};
    }
    self.postMessage({id:data.id, values:call(wasm, data.request).values});
  } catch (error) {
    self.postMessage({id:data.id, error:error.message});
  }
}
// Serialize requests so an asynchronous file read cannot replace another request's cache.
let queue = Promise.resolve();
self.onmessage = ({data}) => { queue = queue.then(() => decode(data)); };
