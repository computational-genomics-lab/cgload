// Node's own fetch, in place of node-fetch v2, which cannot read the data: URL
// the BGZF decoder loads its WebAssembly from.
module.exports = globalThis.fetch;
module.exports.default = globalThis.fetch;
module.exports.Headers = globalThis.Headers;
module.exports.Request = globalThis.Request;
module.exports.Response = globalThis.Response;
