const CACHE = 'scanscope-beta';
const ASSETS = [
  './index.html',
  './apriltag.js',
  './apriltag_wasm.js',
  './apriltag_wasm.wasm',
  './icon.svg',
  './manifest.json',
];

self.addEventListener('install', e => {
  e.waitUntil(caches.open(CACHE).then(c => c.addAll(ASSETS)).then(() => self.skipWaiting()));
});

self.addEventListener('activate', e => {
  e.waitUntil(self.clients.claim());
});

self.addEventListener('fetch', e => {
  const url = new URL(e.request.url);
  if (url.pathname.endsWith('version.json')) {
    e.respondWith(fetch(e.request).catch(() => caches.match(e.request)));
    return;
  }
  e.respondWith(caches.match(e.request).then(r => r || fetch(e.request)));
});

self.addEventListener('message', e => {
  if (e.data && e.data.type === 'UPDATE_NOW') {
    caches.open(CACHE).then(async c => {
      await Promise.all(ASSETS.map(a => fetch(a + '?u=' + Date.now()).then(r => c.put(a, r)).catch(() => {})));
      e.source.postMessage({ type: 'UPDATE_DONE' });
    });
  }
});
