// ④ 視点バイアス解析 — photo_align「🎯視点バイアステスト」のログ末尾 [解析用] JSON を DATA に貼って `node viewbias_test.js`
// 同じ静止ボスを複数視点で撮った採用行（ゲート通過済み）の world yaw が、視点幾何（方位az/画面オフセットoff/距離dist/高さh）で
// どれだけ動くか＝視点バイアスを1変数回帰で炙り出す。傾きが補正係数候補（photo版VIEWBIAS本校正データ）。
//
// 各行: { yaw:世界yaw°, az:ボス基準カメラ方位°, off:ボス中心の画面オフセットpx, dist:距離mm, h:カメラ高さcm }
// ※ボスは静止なので真のyawは一定。yawの変動はすべて測定系の視点バイアス＋ノイズ。

const DATA = [
  // 例: { yaw: 73.08, az: 12, off: 240, dist: 700, h: 113 },
];

function reg(xs, ys){ // 単回帰 y=a+b x、戻り{a,b,r2,rmse}
  const n=xs.length, mx=xs.reduce((s,v)=>s+v,0)/n, my=ys.reduce((s,v)=>s+v,0)/n;
  let sxx=0,sxy=0,syy=0;
  for(let i=0;i<n;i++){ const dx=xs[i]-mx, dy=ys[i]-my; sxx+=dx*dx; sxy+=dx*dy; syy+=dy*dy; }
  const b = sxx>1e-9 ? sxy/sxx : 0, a = my - b*mx;
  let sse=0; for(let i=0;i<n;i++){ const e=ys[i]-(a+b*xs[i]); sse+=e*e; }
  const r2 = syy>1e-9 ? 1-sse/syy : 0;
  return { a, b, r2, rmse: Math.sqrt(sse/n) };
}

if(DATA.length < 3){ console.log('DATA に採用行を3つ以上貼ってください（photo_alignのログ [解析用] JSON）。'); process.exit(0); }

const yaw = DATA.map(d=>d.yaw);
const n = yaw.length, mean = yaw.reduce((s,v)=>s+v,0)/n;
const sd = Math.sqrt(yaw.reduce((s,v)=>s+(v-mean)**2,0)/n);
const spread = Math.max(...yaw)-Math.min(...yaw);

console.log(`採用 ${n}枚`);
console.log(`world yaw: 平均 ${mean.toFixed(3)}°  ばらつき(max-min) ${spread.toFixed(3)}°  σ ${sd.toFixed(3)}°  ＝視点バイアス残差(生)`);
console.log('');

const vars = [
  ['方位 az [°]',  'az',   '°/°',   '°'],
  ['画面offset [px]', 'off', '°/px',  'px'],
  ['距離 dist [mm]',  'dist','°/mm',  'mm'],
  ['高さ h [cm]',     'h',   '°/cm',  'cm'],
];
console.log('── yaw を各視点変数で1変数回帰（傾き=その変数1単位あたりのyaw変化＝バイアス感度）──');
let best=null;
for(const [label,key,unit] of vars){
  const xs=DATA.map(d=>d[key]);
  const r=reg(xs,yaw);
  console.log(`  ${label.padEnd(16)} 傾き ${r.b>=0?'+':''}${r.b.toFixed(4)} ${unit}  R² ${r.r2.toFixed(3)}  残差σ ${r.rmse.toFixed(3)}°`);
  if(best==null || r.r2>best.r2) best={label,key,...r};
}
console.log('');
console.log(`最も効いてる視点変数 = ${best.label}（R² ${best.r2.toFixed(3)}・傾き ${best.b.toFixed(4)}）`);
console.log(`→ これで補正すると残差σ ${best.rmse.toFixed(3)}°（補正前 ${sd.toFixed(3)}°）`);
console.log('');
console.log('判定の目安: 補正後残差σ ≤0.2°なら個別トー実用圏。az(方位)が主因ならphoto版VIEWBIAS_C=その傾きを移植候補。');
console.log('R²が全変数で低い＝単一視点変数で説明できない＝ノイズ支配（ゲートを締める/枚数を増やす/構図を揃える）。');
