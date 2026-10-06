let wasmPromise;
async function module() {
  if (!wasmPromise) wasmPromise = (async () => {
    const response = await fetch('./expresso_viewer.wasm');
    if (!response.ok) throw new Error('WASM decoder missing. Follow the viewer build instructions in the README.');
    return (await WebAssembly.instantiate(await response.arrayBuffer(), {})).instance.exports;
  })();
  return wasmPromise;
}
self.onmessage = async ({data}) => {
  let wasm, requestPtr, vectorPtr, resultPtr, resultLen;
  const request = new TextEncoder().encode(JSON.stringify(data.request));
  try {
    wasm = await module();
    const vector = new Uint8Array(await data.file.slice(data.offset, data.offset + data.length).arrayBuffer());
    if (vector.length !== data.length) throw new Error('Truncated dataset block');
    requestPtr = wasm.allocate(request.length);
    vectorPtr = wasm.allocate(vector.length);
    new Uint8Array(wasm.memory.buffer, requestPtr, request.length).set(request);
    new Uint8Array(wasm.memory.buffer, vectorPtr, vector.length).set(vector);
    const packed = wasm.decode_vector(requestPtr, request.length, vectorPtr, vector.length);
    resultPtr = Number(packed & 0xffffffffn);
    resultLen = Number(packed >> 32n);
    const result = JSON.parse(new TextDecoder().decode(new Uint8Array(wasm.memory.buffer, resultPtr, resultLen)));
    if (result.error) throw new Error(result.error);
    self.postMessage({id:data.id, values:result.values});
  } catch (error) {
    self.postMessage({id:data.id, error:error.message});
  } finally {
    if (wasm) {
      if (requestPtr !== undefined) wasm.release(requestPtr, request.length);
      if (vectorPtr !== undefined) wasm.release(vectorPtr, data.length);
      if (resultPtr !== undefined) wasm.release(resultPtr, resultLen);
    }
  }
};
