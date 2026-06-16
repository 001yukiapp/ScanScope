// 床自己位置: コーナーPnP vs 中心Horn の合成テスト
// 仮説: カメラyaw誤差 ≈ 床タグ中心ノイズ ÷ 床ベースライン。中心Hornは短ベースラインで悪化。
//       コーナーPnPは各タグの四隅(向き)が回転を直接拘束するので短ベースラインでも崩れにくい。
// 実機データ(広い床σ0.9°/寄せた床σ2.6° ≒ ベースライン比2.5倍と一致)の根治策の妥当性を確認する。
const W = require('./phase0/wm_geom.js');
const { hornPose, pnpRefine, projectK, matVec3 } = W;

const FX = 2868, IMGW = 4032, IMGH = 3024;
const K = { fx:FX, fy:FX, cx:IMGW/2, cy:IMGH/2, k1:0 };
const H = 0.105, h = H/2;                 // 小タイル黒部105mm
const LOCAL = [[-h,-h,0],[h,-h,0],[h,h,0],[-h,h,0]];

function mul(A,B){const C=[[0,0,0],[0,0,0],[0,0,0]];for(let i=0;i<3;i++)for(let j=0;j<3;j++)for(let k=0;k<3;k++)C[i][j]+=A[i][k]*B[k][j];return C;}
function T(R){return [[R[0][0],R[1][0],R[2][0]],[R[0][1],R[1][1],R[2][1]],[R[0][2],R[1][2],R[2][2]]];}
function mv(R,v){return [R[0][0]*v[0]+R[0][1]*v[1]+R[0][2]*v[2],R[1][0]*v[0]+R[1][1]*v[1]+R[1][2]*v[2],R[2][0]*v[0]+R[2][1]*v[1]+R[2][2]*v[2]];}
function Rx(a){const c=Math.cos(a),s=Math.sin(a);return[[1,0,0],[0,c,-s],[0,s,c]];}
function Ry(a){const c=Math.cos(a),s=Math.sin(a);return[[c,0,s],[0,1,0],[-s,0,c]];}
// 床タグ(面を上に向けて水平に置く)の world回転: local z(法線)→world +y, in-plane yaw θ
function tagR(th){const c=Math.cos(th),s=Math.sin(th);return [[c,-s,0],[0,0,1],[s,c,0]];}
function geodesic(A,B){const E=mul(A,T(B));const tr=E[0][0]+E[1][1]+E[2][2];return Math.acos(Math.max(-1,Math.min(1,(tr-1)/2)))*180/Math.PI;}
function randn(){let u=0,v=0;while(!u)u=Math.random();while(!v)v=Math.random();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}

// 床タグ配置(world・y=0平面)。spread=中心の散らばり半径
function makeFloor(spread){
  const tags=[];
  const pos=[[-1,-1],[1,-1],[1,1],[-1,1]];  // 4枚を正方
  for(let i=0;i<4;i++){ const th=(i*37)*Math.PI/180; // 各タグ別yaw
    tags.push({ id:44+i, c:[pos[i][0]*spread,0,pos[i][1]*spread], R:tagR(th) }); }
  return tags;
}
function cross(a,b){return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];}
function norm(a){const n=Math.hypot(a[0],a[1],a[2]);return [a[0]/n,a[1]/n,a[2]/n];}
// 真のカメラ world→cam: camPosからtargetを見る(look-at)。cam +z=視線方向
function camPose(camPos, target){
  const f=norm([target[0]-camPos[0],target[1]-camPos[1],target[2]-camPos[2]]);
  const r=norm(cross([0,0,1],f)), u=cross(f,r);
  const Rc2w=[[r[0],u[0],f[0]],[r[1],u[1],f[1]],[r[2],u[2],f[2]]]; // cols=cam軸(world)
  const R=T(Rc2w); const t=mv(R,[-camPos[0],-camPos[1],-camPos[2]]);
  return {R,t,camPos};
}

function trial(spread, sigPx){
  const tags=makeFloor(spread); const cam=camPose([0.15,1.1,-0.35],[0,0,0]);
  // 各タグ: world corners, 投影uv(+ノイズ), 単タグposeから中心obs
  const worldR={}, worldC={}, vis=[];
  for(const tg of tags){
    const Xw = LOCAL.map(l=>{const v=mv(tg.R,l);return [tg.c[0]+v[0],tg.c[1]+v[1],tg.c[2]+v[2]];});
    const uv = Xw.map(X=>{const Y=mv(cam.R,X);const p=projectK(K,[Y[0]+cam.t[0],Y[1]+cam.t[1],Y[2]+cam.t[2]]);return [p.u+randn()*sigPx, p.v+randn()*sigPx];});
    // 単タグposeで中心obs(cam框)を得る(現実の中心ノイズ=深さ増幅を再現)
    const Rtc_true = mul(cam.R, tg.R);              // tag-local→cam 真
    const ttc_true = (()=>{const Y=mv(cam.R,tg.c);return [Y[0]+cam.t[0],Y[1]+cam.t[1],Y[2]+cam.t[2]];})();
    const st = pnpRefine(LOCAL, uv, K, Rtc_true, ttc_true, 12);
    const centerCam = st ? st.t : ttc_true;        // = タグ中心のcam座標(=floorK1T相当)
    worldR[tg.id]=tg.R; worldC[tg.id]=tg.c;
    vis.push({ id:tg.id, uv, centerCam, Xw });
  }
  // --- 中心Horn (現行localizeCamera) ---
  const P=vis.map(v=>worldC[v.id]), C=vis.map(v=>v.centerCam);
  const solH=hornPose(P,C);                          // world→cam
  // --- コーナーPnP (新) : 全タグ四隅でpnpRefine, 初期=中心Horn ---
  const t0=(()=>{const pbar=[0,0,0],cbar=[0,0,0];for(const v of vis){for(let k=0;k<3;k++){pbar[k]+=worldC[v.id][k];cbar[k]+=v.centerCam[k];}}for(let k=0;k<3;k++){pbar[k]/=vis.length;cbar[k]/=vis.length;}const Rp=mv(solH.R,pbar);return [cbar[0]-Rp[0],cbar[1]-Rp[1],cbar[2]-Rp[2]];})();
  const X=[],uvs=[];
  for(const v of vis){ for(let k=0;k<4;k++){ X.push(v.Xw[k]); uvs.push(v.uv[k]); } }
  const solP=pnpRefine(X, uvs, K, solH.R, t0, 15);
  return { errHorn:geodesic(solH.R, cam.R), errPnP: solP?geodesic(solP.R, cam.R):NaN };
}

function run(spread, sigPx, N){
  let sh=0,sp=0; for(let i=0;i<N;i++){const r=trial(spread,sigPx); sh+=r.errHorn; sp+=r.errPnP;}
  return { horn:sh/N, pnp:sp/N };
}

console.log('=== 規約チェック(ノイズ0): 両方とも真姿勢を復元するはず ===');
{ const r=trial(0.7,0); console.log(`  中心Horn誤差 ${r.errHorn.toFixed(4)}°  コーナーPnP誤差 ${r.errPnP.toFixed(4)}°  (両方≈0でPASS)`); }

console.log('\n=== ノイズ1px・ベースライン別の平均回転誤差(30回) ===');
for(const sp of [0.30, 0.50, 0.80, 1.20]){
  const r=run(sp, 1.0, 30);
  console.log(`  床散らばり±${(sp*100).toFixed(0)}cm: 中心Horn ${r.horn.toFixed(3)}°  →  コーナーPnP ${r.pnp.toFixed(3)}°  (改善${(r.horn/r.pnp).toFixed(1)}倍)`);
}
console.log('\n期待: ①ノイズ0で両方≈0(規約OK) ②中心Hornは狭いほど悪化(1/baseline) ③コーナーPnPが狭ベースラインで大幅改善');
