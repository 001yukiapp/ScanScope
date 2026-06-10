// ===== WheelMetry 位置ベース幾何コア (WMGeom) =====
// 傾き(tilt)を使わず、タグ中心の3D位置だけで本体姿勢を出す。
// (1) 距離行列から3D配置を復元（トリラテレーション）
// (2) 復元配置と観測点を剛体整合（Horn法・absolute orientation）
// ブラウザでは window.WMGeom、Nodeでは module.exports に同じAPIを公開。
(function(root){

// --- 3Dベクトル ---
const vsub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]];
const vdot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const vscale=(a,s)=>[a[0]*s,a[1]*s,a[2]*s];
const dvec=(a,b)=>Math.hypot(a[0]-b[0],a[1]-b[1],a[2]-b[2]);

// --- 3x3 線形ソルバ（ガウス消去・部分ピボット）---
function solve3(A,b){
    const M=[[A[0][0],A[0][1],A[0][2],b[0]],
             [A[1][0],A[1][1],A[1][2],b[1]],
             [A[2][0],A[2][1],A[2][2],b[2]]];
    for(let c=0;c<3;c++){
        let p=c; for(let r=c+1;r<3;r++) if(Math.abs(M[r][c])>Math.abs(M[p][c])) p=r;
        if(Math.abs(M[p][c])<1e-12) return null;
        [M[c],M[p]]=[M[p],M[c]];
        for(let r=0;r<3;r++){ if(r===c) continue; const f=M[r][c]/M[c][c];
            for(let k=c;k<4;k++) M[r][k]-=f*M[c][k]; }
    }
    return [M[0][3]/M[0][0], M[1][3]/M[1][1], M[2][3]/M[2][2]];
}

// --- 距離行列から3D配置を復元 ---
// ids: 配列, dist(i,j): 既知距離 or null/0
// 返り値 { pos:{id:[x,y,z]}, placed:[id], failed:[id] }
function reconstruct(ids, dist){
    const pos={};
    // 1) 非退化な種三角形を探す
    let seed=null;
    outer:
    for(let a=0;a<ids.length;a++)for(let b=a+1;b<ids.length;b++)for(let c=b+1;c<ids.length;c++){
        const i=ids[a],j=ids[b],k=ids[c];
        const dij=dist(i,j),dik=dist(i,k),djk=dist(j,k);
        if(dij&&dik&&djk){
            const x=(dij*dij+dik*dik-djk*djk)/(2*dij);
            const y2=dik*dik-x*x;
            if(y2>1e-6){ seed={i,j,k,dij,x,y:Math.sqrt(y2)}; break outer; }
        }
    }
    if(!seed) return { pos:{}, placed:[], failed:ids.slice() };
    pos[seed.i]=[0,0,0]; pos[seed.j]=[seed.dij,0,0]; pos[seed.k]=[seed.x,seed.y,0];
    const placed=new Set([seed.i,seed.j,seed.k]);

    // 2) 4点目: 種3点(同一平面z=0)への距離から ±z を二次式で解き +z を採用（カイラリティ固定）
    function placeQuadratic(k){
        const anc=[...placed].filter(p=>dist(k,p)).slice(0,3);
        if(anc.length<3) return false;
        const P0=pos[anc[0]], d0=dist(k,anc[0]);
        const rows=[],rhs=[];
        for(let t=1;t<3;t++){ const Pa=pos[anc[t]], da=dist(k,anc[t]);
            rows.push([2*(Pa[0]-P0[0]),2*(Pa[1]-P0[1])]);
            rhs.push((vdot(Pa,Pa)-vdot(P0,P0))-(da*da-d0*d0));
        }
        const det=rows[0][0]*rows[1][1]-rows[0][1]*rows[1][0];
        if(Math.abs(det)<1e-9) return false;
        const x=(rhs[0]*rows[1][1]-rhs[1]*rows[0][1])/det;
        const y=(rows[0][0]*rhs[1]-rows[1][0]*rhs[0])/det;
        const z2=d0*d0-(x-P0[0])**2-(y-P0[1])**2-(0-P0[2])**2;
        const z=z2>0?Math.sqrt(z2):0;
        pos[k]=[x,y,z]; placed.add(k); return true;
    }
    for(const k of ids){ if(placed.has(k)) continue; if(placeQuadratic(k)) break; }

    // 3) 残り: 非同一平面の配置点≥4から線形最小二乗で一意決定
    function placeLinear(k){
        const anc=[...placed].filter(p=>dist(k,p));
        if(anc.length<4) return false;
        const P0=pos[anc[0]], d0=dist(k,anc[0]);
        const ATA=[[0,0,0],[0,0,0],[0,0,0]], ATb=[0,0,0];
        for(let t=1;t<anc.length;t++){ const Pa=pos[anc[t]], da=dist(k,anc[t]);
            const row=[2*(Pa[0]-P0[0]),2*(Pa[1]-P0[1]),2*(Pa[2]-P0[2])];
            const r=(vdot(Pa,Pa)-vdot(P0,P0))-(da*da-d0*d0);
            for(let a=0;a<3;a++){ for(let b=0;b<3;b++) ATA[a][b]+=row[a]*row[b]; ATb[a]+=row[a]*r; }
        }
        const sol=solve3(ATA,ATb);
        if(!sol) return false;
        pos[k]=sol; placed.add(k); return true;
    }
    let progress=true;
    while(progress){ progress=false;
        for(const k of ids){ if(placed.has(k)) continue; if(placeLinear(k)){ progress=true; } }
    }
    const failed=ids.filter(k=>!placed.has(k));
    return { pos, placed:[...placed], failed };
}

// --- 4x4対称行列のヤコビ固有値分解 ---
function jacobiEigen4(Ain){
    const n=4, A=Ain.map(r=>r.slice());
    const V=[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
    for(let sweep=0;sweep<100;sweep++){
        let off=0; for(let p=0;p<n;p++)for(let q=p+1;q<n;q++) off+=A[p][q]*A[p][q];
        if(off<1e-20) break;
        for(let p=0;p<n;p++)for(let q=p+1;q<n;q++){
            if(Math.abs(A[p][q])<1e-18) continue;
            const phi=0.5*Math.atan2(2*A[p][q], A[q][q]-A[p][p]);
            const c=Math.cos(phi), s=Math.sin(phi);
            for(let k=0;k<n;k++){ const akp=A[k][p], akq=A[k][q]; A[k][p]=c*akp-s*akq; A[k][q]=s*akp+c*akq; }
            for(let k=0;k<n;k++){ const apk=A[p][k], aqk=A[q][k]; A[p][k]=c*apk-s*aqk; A[q][k]=s*apk+c*aqk; }
            for(let k=0;k<n;k++){ const vkp=V[k][p], vkq=V[k][q]; V[k][p]=c*vkp-s*vkq; V[k][q]=s*vkp+c*vkq; }
        }
    }
    return { vals:[A[0][0],A[1][1],A[2][2],A[3][3]], vecs:V };
}

const quatToMat=q=>{ const [w,x,y,z]=q; return [
    [1-2*(y*y+z*z), 2*(x*y-w*z),   2*(x*z+w*y)],
    [2*(x*y+w*z),   1-2*(x*x+z*z), 2*(y*z-w*x)],
    [2*(x*z-w*y),   2*(y*z+w*x),   1-2*(x*x+y*y)] ]; };

// --- Horn法: 対応点 P(本体) ↔ C(カメラ) から回転 R(本体→カメラ) を解く ---
// w: オプション重み配列（省略時は等重み）
function hornPose(P, C, w){
    const n=P.length;
    let wsum=0; const wt=[];
    for(let i=0;i<n;i++){ wt.push(w?w[i]:1); wsum+=wt[i]; }
    let pbar=[0,0,0], cbar=[0,0,0];
    for(let i=0;i<n;i++){ for(let a=0;a<3;a++){ pbar[a]+=wt[i]*P[i][a]; cbar[a]+=wt[i]*C[i][a]; } }
    pbar=vscale(pbar,1/wsum); cbar=vscale(cbar,1/wsum);
    const S=[[0,0,0],[0,0,0],[0,0,0]];
    for(let i=0;i<n;i++){ const p=vsub(P[i],pbar), c=vsub(C[i],cbar);
        for(let a=0;a<3;a++)for(let b=0;b<3;b++) S[a][b]+=wt[i]*p[a]*c[b]; }
    const Sxx=S[0][0],Sxy=S[0][1],Sxz=S[0][2],
          Syx=S[1][0],Syy=S[1][1],Syz=S[1][2],
          Szx=S[2][0],Szy=S[2][1],Szz=S[2][2];
    const N=[
        [Sxx+Syy+Szz, Syz-Szy,      Szx-Sxz,      Sxy-Syx     ],
        [Syz-Szy,     Sxx-Syy-Szz,  Sxy+Syx,      Szx+Sxz     ],
        [Szx-Sxz,     Sxy+Syx,     -Sxx+Syy-Szz,  Syz+Szy     ],
        [Sxy-Syx,     Szx+Sxz,      Syz+Szy,     -Sxx-Syy+Szz ],
    ];
    const { vals, vecs }=jacobiEigen4(N);
    let mi=0; for(let i=1;i<4;i++) if(vals[i]>vals[mi]) mi=i;
    let q=[vecs[0][mi],vecs[1][mi],vecs[2][mi],vecs[3][mi]];
    const qn=Math.hypot(q[0],q[1],q[2],q[3])||1; q=q.map(v=>v/qn);
    const R=quatToMat(q);
    let se=0; for(let i=0;i<n;i++){ const p=vsub(P[i],pbar), c=vsub(C[i],cbar);
        const Rp=[R[0][0]*p[0]+R[0][1]*p[1]+R[0][2]*p[2],
                  R[1][0]*p[0]+R[1][1]*p[1]+R[1][2]*p[2],
                  R[2][0]*p[0]+R[2][1]*p[1]+R[2][2]*p[2]];
        se+=wt[i]*((Rp[0]-c[0])**2+(Rp[1]-c[1])**2+(Rp[2]-c[2])**2); }
    return { R, q, rms:Math.sqrt(se/wsum) };
}

// --- 外れ値除去付きHorn（残差median+2.5*MAD超えを除外して再推定）---
// w: オプション重み配列（入射角cosなど）
function hornPoseRobust(P, C, w){
    if(P.length<3) return {...hornPose(P,C,w), outliers:0};
    const sol0=hornPose(P,C,w);
    // 重み付き重心でper-point残差を計算
    const wt=P.map((_,i)=>w?w[i]:1), wsum=wt.reduce((a,b)=>a+b,0);
    let pbar=[0,0,0], cbar=[0,0,0];
    for(let i=0;i<P.length;i++) for(let a=0;a<3;a++){ pbar[a]+=wt[i]*P[i][a]; cbar[a]+=wt[i]*C[i][a]; }
    pbar=vscale(pbar,1/wsum); cbar=vscale(cbar,1/wsum);
    const R=sol0.R;
    const res=P.map((Pi,i)=>{
        const p=vsub(Pi,pbar), c=vsub(C[i],cbar);
        const Rp=[R[0][0]*p[0]+R[0][1]*p[1]+R[0][2]*p[2],
                  R[1][0]*p[0]+R[1][1]*p[1]+R[1][2]*p[2],
                  R[2][0]*p[0]+R[2][1]*p[1]+R[2][2]*p[2]];
        return Math.sqrt((Rp[0]-c[0])**2+(Rp[1]-c[1])**2+(Rp[2]-c[2])**2);
    });
    const srt=[...res].sort((a,b)=>a-b);
    const med=srt[Math.floor(srt.length/2)];
    const mad=([...res.map(r=>Math.abs(r-med))].sort((a,b)=>a-b))[Math.floor(res.length/2)]*1.4826;
    const th=med+2.5*Math.max(mad,5e-4);
    const mask=res.map(r=>r<=th);
    const nOut=mask.filter(v=>!v).length;
    const nKeep=P.length-nOut;
    if(nOut===0||nKeep<3) return {...sol0, outliers:0};
    const Pk=P.filter((_,i)=>mask[i]), Ck=C.filter((_,i)=>mask[i]);
    const wk=w?w.filter((_,i)=>mask[i]):null;
    return {...hornPose(Pk,Ck,wk), outliers:nOut};
}

// ===== コーナー再投影PnP（タグ4隅で姿勢を直接最適化）=====
// 中心位置ベースHornの弱点（タグ中心の奥行きバイアスが小スパン本体で回転に化ける）への対策。
// 画像コーナー座標を直接使うため奥行きバイアスの影響を受けない。

const matT3=R=>[[R[0][0],R[1][0],R[2][0]],[R[0][1],R[1][1],R[2][1]],[R[0][2],R[1][2],R[2][2]]];
function matMul3(A,B){
    const C=[[0,0,0],[0,0,0],[0,0,0]];
    for(let i=0;i<3;i++)for(let j=0;j<3;j++){let s=0;for(let k=0;k<3;k++)s+=A[i][k]*B[k][j];C[i][j]=s;}
    return C;
}
const matVec3=(R,v)=>[
    R[0][0]*v[0]+R[0][1]*v[1]+R[0][2]*v[2],
    R[1][0]*v[0]+R[1][1]*v[1]+R[1][2]*v[2],
    R[2][0]*v[0]+R[2][1]*v[1]+R[2][2]*v[2]];

// --- 汎用NxN線形ソルバ（ガウス消去・部分ピボット）---
function solveN(A,b){
    const n=b.length, M=A.map((r,i)=>[...r,b[i]]);
    for(let c=0;c<n;c++){
        let p=c; for(let r=c+1;r<n;r++) if(Math.abs(M[r][c])>Math.abs(M[p][c])) p=r;
        if(Math.abs(M[p][c])<1e-12) return null;
        [M[c],M[p]]=[M[p],M[c]];
        for(let r=0;r<n;r++){ if(r===c) continue; const f=M[r][c]/M[c][c];
            for(let k=c;k<=n;k++) M[r][k]-=f*M[c][k]; }
    }
    const x=[]; for(let i=0;i<n;i++) x.push(M[i][n]/M[i][i]);
    return x;
}

// --- SO(3)指数写像（Rodrigues）---
function expSO3(w){
    const th=Math.hypot(w[0],w[1],w[2]);
    if(th<1e-12) return [[1,0,0],[0,1,0],[0,0,1]];
    const k=[w[0]/th,w[1]/th,w[2]/th], c=Math.cos(th), s=Math.sin(th), C=1-c;
    return [
        [c+k[0]*k[0]*C,      k[0]*k[1]*C-k[2]*s, k[0]*k[2]*C+k[1]*s],
        [k[1]*k[0]*C+k[2]*s, c+k[1]*k[1]*C,      k[1]*k[2]*C-k[0]*s],
        [k[2]*k[0]*C-k[1]*s, k[2]*k[1]*C+k[0]*s, c+k[2]*k[2]*C]];
}

// --- 回転間の角度（度）---
function rotAngleDeg(A,B){
    const R=matMul3(matT3(A),B);
    const tr=R[0][0]+R[1][1]+R[2][2];
    return Math.acos(Math.max(-1,Math.min(1,(tr-1)/2)))*180/Math.PI;
}

// --- 近接回転の平均（行列成分平均→列Gram-Schmidt直交化）---
function avgRotGS(Rs){
    const M=[[0,0,0],[0,0,0],[0,0,0]];
    for(const R of Rs) for(let i=0;i<3;i++)for(let j=0;j<3;j++) M[i][j]+=R[i][j];
    const nrm=v=>{const n=Math.hypot(v[0],v[1],v[2]); return n>1e-12?[v[0]/n,v[1]/n,v[2]/n]:v;};
    let c0=nrm([M[0][0],M[1][0],M[2][0]]);
    let c1=[M[0][1],M[1][1],M[2][1]];
    const d=c0[0]*c1[0]+c0[1]*c1[1]+c0[2]*c1[2];
    c1=nrm([c1[0]-c0[0]*d, c1[1]-c0[1]*d, c1[2]-c0[2]*d]);
    const c2=[c0[1]*c1[2]-c0[2]*c1[1], c0[2]*c1[0]-c0[0]*c1[2], c0[0]*c1[1]-c0[1]*c1[0]];
    return [[c0[0],c1[0],c2[0]],[c0[1],c1[1],c2[1]],[c0[2],c1[2],c2[2]]];
}

// --- 候補回転群から最大クラスタの代表を選ぶ（flip解の排除）---
// cands: [R,...]  weights: 任意の事前重み（凸形状の法線外向き度など。省略時等重み）
// 戻り { R, support(0〜1), spreadDeg } | null
function pickStableRotation(cands, thDeg, weights){
    if(!cands||!cands.length) return null;
    const th=thDeg||15;
    const w=weights||cands.map(()=>1);
    const totalW=w.reduce((a,b)=>a+b,0)||1;
    let best=null, bestW=-1;
    for(let i=0;i<cands.length;i++){
        const idx=[];
        for(let j=0;j<cands.length;j++) if(rotAngleDeg(cands[i],cands[j])<th) idx.push(j);
        const cw=idx.reduce((s,j)=>s+w[j],0);
        if(cw>bestW){ bestW=cw; best=idx; }
    }
    const R=avgRotGS(best.map(j=>cands[j]));
    let spread=0; for(const j of best) spread+=rotAngleDeg(R,cands[j]);
    return { R, support:bestW/totalW, spreadDeg:spread/best.length };
}

// --- コーナー再投影PnP（Gauss-Newton・6DoF）---
// X: [[x,y,z]…] 本体系3D点 / uv: [[u,v]…] 観測px / K:{fx,fy,cx,cy}
// R0,t0: 初期姿勢（Horn等から）。戻り { R, t, rmsPx, iters } | null
function pnpRefine(X, uv, K, R0, t0, maxIters){
    let R=R0.map(r=>r.slice()), t=t0.slice();
    const iters=maxIters||10;
    let rms=Infinity, it=0;
    for(it=0; it<iters; it++){
        const H=Array.from({length:6},()=>new Array(6).fill(0));
        const g=new Array(6).fill(0);
        let se=0;
        for(let i=0;i<X.length;i++){
            const Z=matVec3(R,X[i]);
            const Y=[Z[0]+t[0], Z[1]+t[1], Z[2]+t[2]];
            if(Y[2]<1e-6) return null;
            const iz=1/Y[2];
            const ru=K.fx*Y[0]*iz+K.cx-uv[i][0];
            const rv=K.fy*Y[1]*iz+K.cy-uv[i][1];
            se+=ru*ru+rv*rv;
            const du=[K.fx*iz, 0, -K.fx*Y[0]*iz*iz];
            const dv=[0, K.fy*iz, -K.fy*Y[1]*iz*iz];
            // dY/dδw = -skew(Z)（左摂動 R←exp(δw)R）, dY/dδt = I
            const Jw=[[0,Z[2],-Z[1]],[-Z[2],0,Z[0]],[Z[1],-Z[0],0]];
            const Ju=[du[0]*Jw[0][0]+du[1]*Jw[1][0]+du[2]*Jw[2][0],
                      du[0]*Jw[0][1]+du[1]*Jw[1][1]+du[2]*Jw[2][1],
                      du[0]*Jw[0][2]+du[1]*Jw[1][2]+du[2]*Jw[2][2],
                      du[0],du[1],du[2]];
            const Jv=[dv[0]*Jw[0][0]+dv[1]*Jw[1][0]+dv[2]*Jw[2][0],
                      dv[0]*Jw[0][1]+dv[1]*Jw[1][1]+dv[2]*Jw[2][1],
                      dv[0]*Jw[0][2]+dv[1]*Jw[1][2]+dv[2]*Jw[2][2],
                      dv[0],dv[1],dv[2]];
            for(let a=0;a<6;a++){
                g[a]+=Ju[a]*ru+Jv[a]*rv;
                for(let b=a;b<6;b++) H[a][b]+=Ju[a]*Ju[b]+Jv[a]*Jv[b];
            }
        }
        for(let a=0;a<6;a++)for(let b=0;b<a;b++) H[a][b]=H[b][a];
        for(let a=0;a<6;a++) H[a][a]+=1e-9;          // 数値安定化
        const newRms=Math.sqrt(se/X.length);
        const d=solveN(H, g.map(x=>-x));
        if(!d) break;
        R=matMul3(expSO3([d[0],d[1],d[2]]), R);
        t=[t[0]+d[3], t[1]+d[4], t[2]+d[5]];
        if(Math.abs(rms-newRms)<1e-4){ rms=newRms; break; }
        rms=newRms;
    }
    return { R, t, rmsPx:rms, iters:it };
}

// --- 鏡像（z反転）---
function mirrorPos(pos){ const o={}; for(const k in pos){ const p=pos[k]; o[k]=[p[0],p[1],-p[2]]; } return o; }
// --- 配置が観測フレーム群にどれだけ合うか（平均RMS）---
function chiralityResidual(pos, frames){
    let s=0,n=0;
    for(const fr of frames){
        const ids=Object.keys(fr).filter(id=>pos[id]);
        if(ids.length<4) continue;
        const P=ids.map(id=>pos[id]), C=ids.map(id=>fr[id]);
        s+=hornPose(P,C).rms; n++;
    }
    return n? s/n : Infinity;
}
// --- キャリブ観測フレーム群から本体配置を構築（距離復元＋カイラリティ確定）---
// frames: [ {id:[x,y,z](カメラ系・メートル), ...}, ... ]
function buildBody(ids, frames){
    const obs={}, key=(i,j)=>i<j?i+','+j:j+','+i;
    for(const fr of frames){ const fids=Object.keys(fr).map(Number);
        for(let a=0;a<fids.length;a++)for(let b=a+1;b<fids.length;b++){
            const i=fids[a],j=fids[b]; (obs[key(i,j)]=obs[key(i,j)]||[]).push(dvec(fr[i],fr[j]));
        } }
    const med=arr=>{ const s=[...arr].sort((x,y)=>x-y); return s[Math.floor(s.length/2)]; };
    const dist=(i,j)=> i===j?0:(obs[key(i,j)]?med(obs[key(i,j)]):null);
    let { pos, placed, failed }=reconstruct(ids, dist);
    const m=mirrorPos(pos);
    if(chiralityResidual(m,frames) < chiralityResidual(pos,frames)) pos=m;
    return { pos, placed, failed };
}

const api={ reconstruct, hornPose, hornPoseRobust, buildBody, chiralityResidual, mirrorPos, jacobiEigen4, quatToMat, solve3, dvec,
    pnpRefine, expSO3, rotAngleDeg, avgRotGS, pickStableRotation, solveN, matMul3, matT3, matVec3 };
if(typeof module!=='undefined' && module.exports) module.exports=api;
else root.WMGeom=api;

})(typeof self!=='undefined'?self:this);
