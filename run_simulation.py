#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Monte Carlo Raytracing 小行星粒子模擬 — 完整回測 + PDE 生成
產出：7 張分析圖 + 控制台文字說明
"""

import numpy as np
from scipy.sparse import lil_matrix
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os, time

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ═══════════════════════════════════════════
# STEP 0: Parameters
# ═══════════════════════════════════════════
np.random.seed(42)
SUBDIV      = 3          # level 3 for speed (642 verts, 1280 faces)
TRIAXIAL    = [1.0, 0.8, 0.6]
N_CRATERS   = 30
N_LUMPS     = 8
ROUGHNESS   = 0.005
TT          = 24         # orbital time steps
NMONTE      = 30         # MC samples per face
OBLIQUITY   = 45.0       # degrees

print("=" * 60)
print("  Monte Carlo Raytracing — 小行星粒子模擬回測")
print("=" * 60)

# ═══════════════════════════════════════════
# STEP 1: Icosphere Generation
# ═══════════════════════════════════════════
print("\n[Step 1] 生成 Icosphere 基礎網格...")

def make_icosphere(level):
    t = (1 + np.sqrt(5)) / 2
    V = np.array([
        [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
        [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
        [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]
    ], dtype=float)
    V /= np.linalg.norm(V, axis=1, keepdims=True)

    F = np.array([
        [0,11,5],[0,5,1],[0,1,7],[0,7,10],[0,10,11],
        [1,5,9],[5,11,4],[11,10,2],[10,7,6],[7,1,8],
        [3,9,4],[3,4,2],[3,2,6],[3,6,8],[3,8,9],
        [4,9,5],[2,4,11],[6,2,10],[8,6,7],[9,8,1]
    ], dtype=int)

    for _ in range(level):
        edge_map = {}
        nV = len(V)
        new_F = []
        new_verts = list(V)

        def get_mid(i, j):
            key = (min(i,j), max(i,j))
            if key in edge_map:
                return edge_map[key]
            mid = (new_verts[i] + new_verts[j]) / 2
            mid /= np.linalg.norm(mid)
            idx = len(new_verts)
            new_verts.append(mid)
            edge_map[key] = idx
            return idx

        for f in F:
            a, b, c = f
            ab = get_mid(a, b)
            bc = get_mid(b, c)
            ca = get_mid(c, a)
            new_F.extend([[a,ab,ca],[b,bc,ab],[c,ca,bc],[ab,bc,ca]])

        V = np.array(new_verts)
        F = np.array(new_F, dtype=int)

    return V, F

V, F = make_icosphere(SUBDIV)
print(f"  Icosphere: {len(V)} vertices, {len(F)} faces")

# ═══════════════════════════════════════════
# STEP 2: Triaxial Deformation
# ═══════════════════════════════════════════
print(f"\n[Step 2] 三軸橢球變形 a:b:c = {TRIAXIAL}...")
V_sphere = V.copy()
V[:, 0] *= TRIAXIAL[0]
V[:, 1] *= TRIAXIAL[1]
V[:, 2] *= TRIAXIAL[2]

# ═══════════════════════════════════════════
# STEP 3: Impact Craters
# ═══════════════════════════════════════════
print(f"\n[Step 3] 添加 {N_CRATERS} 個撞擊坑...")
R = np.mean(np.linalg.norm(V, axis=1))

d_min, d_max = 0.05, 0.4
u_max = 1 - (d_min / d_max) ** 2
u = np.random.rand(N_CRATERS) * u_max
diameters = d_min / np.sqrt(1 - u) * R
diameters = np.sort(diameters)[::-1]

centres = np.random.randn(N_CRATERS, 3)
centres /= np.linalg.norm(centres, axis=1, keepdims=True)

V_pre_crater = V.copy()
for ci in range(N_CRATERS):
    D = diameters[ci]
    rad = D / 2
    dep = 0.2 * D
    rim_h = dep * 0.15
    ctr = centres[ci]

    for vi in range(len(V)):
        p = V[vi]
        r_v = np.linalg.norm(p)
        p_hat = p / r_v
        cos_ang = np.clip(np.dot(p_hat, ctr), -1, 1)
        ang = np.arccos(cos_ang)
        dist = ang * R

        if dist < rad * 1.5:
            x = dist / rad
            if x <= 1.0:
                delta = -dep * (1 - x * x)
            elif x <= 1.5:
                t = (x - 1.0) / 0.5
                delta = rim_h * 0.5 * (1 + np.cos(np.pi * t))
            else:
                delta = 0
            V[vi] = p_hat * (r_v + delta)

# ═══════════════════════════════════════════
# STEP 4: Accretion Lumps
# ═══════════════════════════════════════════
print(f"\n[Step 4] 添加 {N_LUMPS} 個增積隆起...")
lump_centres = np.random.randn(N_LUMPS, 3)
lump_centres /= np.linalg.norm(lump_centres, axis=1, keepdims=True)
amps = 0.02 + 0.06 * np.random.rand(N_LUMPS)
sigmas = 0.2 + 0.4 * np.random.rand(N_LUMPS)

for li in range(N_LUMPS):
    ctr = lump_centres[li]
    amp = amps[li] * R
    sigma = sigmas[li]
    for vi in range(len(V)):
        p = V[vi]
        r_v = np.linalg.norm(p)
        p_hat = p / r_v
        ang = np.arccos(np.clip(np.dot(p_hat, ctr), -1, 1))
        delta = amp * np.exp(-ang**2 / (2 * sigma**2))
        if delta > 1e-6 * R:
            V[vi] = p_hat * (r_v + delta)

# ═══════════════════════════════════════════
# STEP 5: Surface Roughness
# ═══════════════════════════════════════════
print(f"\n[Step 5] 添加表面粗糙度 (Laplacian smoothing)...")
nV = len(V)
R = np.mean(np.linalg.norm(V, axis=1))
noise = ROUGHNESS * R * np.random.randn(nV)

# Build adjacency
A = lil_matrix((nV, nV))
for f in F:
    for i in range(3):
        a, b = f[i], f[(i+1)%3]
        A[a, b] = 1
        A[b, a] = 1
A = A.tocsr()

# Laplacian smoothing
for _ in range(3):
    noise_new = np.zeros(nV)
    for vi in range(nV):
        nbrs = A[vi].nonzero()[1]
        if len(nbrs) > 0:
            noise_new[vi] = np.mean(noise[nbrs])
        else:
            noise_new[vi] = noise[vi]
    noise = noise_new

for vi in range(nV):
    r_v = np.linalg.norm(V[vi])
    V[vi] = V[vi] / r_v * (r_v + noise[vi])

# Fix face winding
print("\n[Step 5b] 修正法向量方向...")
n_flipped = 0
for fi in range(len(F)):
    a, b, c = V[F[fi, 0]], V[F[fi, 1]], V[F[fi, 2]]
    centroid = (a + b + c) / 3
    normal = np.cross(c - a, b - a)
    if np.dot(normal, centroid) < 0:
        F[fi] = [F[fi, 0], F[fi, 2], F[fi, 1]]
        n_flipped += 1
print(f"  翻轉了 {n_flipped} 個面")

# Compute face normals and centroids
nface = len(F)
face_normals = np.zeros((nface, 3))
face_centroids = np.zeros((nface, 3))
face_areas = np.zeros(nface)

for fi in range(nface):
    a, b, c = V[F[fi, 0]], V[F[fi, 1]], V[F[fi, 2]]
    face_centroids[fi] = (a + b + c) / 3
    c1 = b - a
    c2 = c - a
    normal = np.cross(c2, c1)
    norm_len = np.linalg.norm(normal)
    if norm_len > 0:
        face_normals[fi] = normal / norm_len
    face_areas[fi] = norm_len / 2

print(f"\n  最終網格: {len(V)} vertices, {nface} faces")
print(f"  平均面積: {np.mean(face_areas):.6f}")

# ═══════════════════════════════════════════
# FIGURE 1: Asteroid Mesh 3D View
# ═══════════════════════════════════════════
print("\n[Figure 1] 繪製小行星 3D 網格...")
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Color by radius
radii = np.linalg.norm(V, axis=1)
face_radii = np.mean(radii[F], axis=1)
face_colors = cm.terrain((face_radii - face_radii.min()) / (face_radii.max() - face_radii.min()))

polys = [V[f] for f in F]
mesh = Poly3DCollection(polys, facecolors=face_colors, edgecolors='k', linewidths=0.1, alpha=0.9)
ax.add_collection3d(mesh)

lim = np.max(np.abs(V)) * 1.2
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
ax.set_title(f'Step 1-5: Asteroid Mesh ({len(V)} verts, {nface} faces)\n'
             f'Triaxial [{TRIAXIAL}], {N_CRATERS} craters, {N_LUMPS} lumps',
             fontsize=12)
ax.view_init(elev=20, azim=30)
fig.savefig(os.path.join(OUT_DIR, 'fig1_asteroid_mesh.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  -> fig1_asteroid_mesh.png saved")

# ═══════════════════════════════════════════
# FIGURE 2: Crater Profile Cross-Section
# ═══════════════════════════════════════════
print("\n[Figure 2] 繪製撞擊坑剖面圖...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Crater profile
x = np.linspace(0, 2, 500)
D_ex = 0.3
dep = 0.2 * D_ex
rim_h = dep * 0.15
profile = np.zeros_like(x)
for i, xi in enumerate(x):
    if xi <= 1.0:
        profile[i] = -dep * (1 - xi**2)
    elif xi <= 1.5:
        t = (xi - 1.0) / 0.5
        profile[i] = rim_h * 0.5 * (1 + np.cos(np.pi * t))

ax1.fill_between(x, profile, -dep*1.2, alpha=0.3, color='brown')
ax1.plot(x, profile, 'k-', linewidth=2)
ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax1.axvline(x=1.0, color='r', linestyle=':', alpha=0.5, label='crater rim (x=1)')
ax1.set_xlabel('Normalized distance (x = dist/radius)')
ax1.set_ylabel('Radial displacement')
ax1.set_title('Crater Profile: Parabolic Bowl + Raised Rim')
ax1.legend()
ax1.annotate(f'depth = {dep:.3f}', xy=(0, -dep), fontsize=10, color='blue')
ax1.annotate(f'rim_h = {rim_h:.4f}', xy=(1.1, rim_h*0.8), fontsize=10, color='red')

# Power-law size distribution
D_vals = np.linspace(0.05, 0.4, 100)
N_gt = (0.05 / D_vals) ** 2 * N_CRATERS
ax2.loglog(D_vals * R, N_gt, 'b-', linewidth=2)
ax2.set_xlabel('Crater Diameter')
ax2.set_ylabel('N(>D) Cumulative Count')
ax2.set_title(r'Crater Size Distribution: $N(>D) \propto D^{-2}$')
ax2.grid(True, alpha=0.3)

fig.suptitle('Step 3: Impact Crater Model', fontsize=14, fontweight='bold')
fig.savefig(os.path.join(OUT_DIR, 'fig2_crater_profile.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  -> fig2_crater_profile.png saved")

# ═══════════════════════════════════════════
# STEP 6: Monte Carlo Shadow + Irradiance
# ═══════════════════════════════════════════
print(f"\n[Step 6] Monte Carlo 光線追蹤 (TT={TT}, Nmonte={NMONTE})...")
za = np.radians(OBLIQUITY)
Flux = np.zeros((nface, TT))
initial_sc = 0.89

t_start = time.time()

for t_idx in range(TT):
    fi = (t_idx + 0.5) / TT * 2 * np.pi
    solar = np.array([np.cos(-za)*np.sin(fi), np.cos(-za)*np.cos(fi), np.sin(-za)])

    for nf in range(nface):
        # Lambert cosine
        c6 = np.dot(solar, face_normals[nf])
        if c6 <= 0:
            Flux[nf, t_idx] = 0
            continue

        # Phase 1: Find shadow neighbors
        shadow_nbrs = []
        s2 = face_centroids[nf]
        for nnf in range(nface):
            if nnf == nf:
                continue
            s1 = face_centroids[nnf]
            line_s = s1 - s2
            line_norm = np.linalg.norm(line_s)
            if line_norm < 1e-10:
                continue
            line_s /= line_norm
            if np.dot(line_s, solar) >= initial_sc:
                shadow_nbrs.append(nnf)

        if len(shadow_nbrs) == 0:
            Flux[nf, t_idx] = c6
            continue

        # Phase 2: Monte Carlo shadow sampling
        a = V[F[nf, 0]]
        b = V[F[nf, 1]]
        c = V[F[nf, 2]]
        occlusion = 0

        for _ in range(NMONTE):
            # Random point on triangle (barycentric)
            r1, r2 = np.random.rand(), np.random.rand()
            if r1 + r2 > 1:
                r1, r2 = 1 - r1, 1 - r2
            pt = a + r1 * (b - a) + r2 * (c - a)

            # Check against each shadow neighbor
            hit = False
            for nnf in shadow_nbrs:
                da = V[F[nnf, 0]]
                db = V[F[nnf, 1]]
                dc = V[F[nnf, 2]]
                fn = face_normals[nnf]

                denom = np.dot(fn, solar)
                if abs(denom) < 1e-12:
                    continue
                t_hit = (np.dot(fn, da) - np.dot(fn, pt)) / denom
                if t_hit <= 0:
                    continue

                endc = pt + solar * t_hit

                # Barycentric area test
                area_full = np.linalg.norm(np.cross(db - da, dc - da))
                a1 = np.linalg.norm(np.cross(endc - da, endc - db))
                a2 = np.linalg.norm(np.cross(endc - da, endc - dc))
                a3 = np.linalg.norm(np.cross(endc - dc, endc - db))
                if (a1 + a2 + a3) - area_full <= 0.001:
                    hit = True
                    break

            if hit:
                occlusion += 1

        Flux[nf, t_idx] = c6 * (1 - occlusion / NMONTE)

    if (t_idx + 1) % 6 == 0:
        elapsed = time.time() - t_start
        print(f"  Time step {t_idx+1}/{TT} done ({elapsed:.1f}s)")

total_time = time.time() - t_start
print(f"  光線追蹤完成！總耗時: {total_time:.1f}s")

# ═══════════════════════════════════════════
# FIGURE 3: Irradiance Map at 4 Time Steps
# ═══════════════════════════════════════════
print("\n[Figure 3] 繪製四個時刻的照度分布...")
fig = plt.figure(figsize=(16, 12))
time_steps = [0, TT//4, TT//2, 3*TT//4]

for idx, t_idx in enumerate(time_steps):
    ax = fig.add_subplot(2, 2, idx+1, projection='3d')
    flux_t = Flux[:, t_idx]

    face_colors = cm.hot(flux_t / max(flux_t.max(), 0.01))
    polys = [V[f] for f in F]
    mesh = Poly3DCollection(polys, facecolors=face_colors, edgecolors='none', alpha=0.95)
    ax.add_collection3d(mesh)

    # Solar direction marker
    fi = (t_idx + 0.5) / TT * 2 * np.pi
    solar_pos = 2.0 * np.array([np.cos(za)*np.sin(fi+np.pi), np.cos(za)*np.cos(fi+np.pi), np.sin(za)])
    ax.scatter(*solar_pos, color='yellow', s=200, marker='*', zorder=10, edgecolors='orange')

    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    phase = t_idx / TT * 360
    ax.set_title(f'Orbital Phase = {phase:.0f} deg\nMax Flux = {flux_t.max():.3f}', fontsize=11)
    ax.view_init(elev=20, azim=30)

fig.suptitle(f'Step 6: Monte Carlo Irradiance (Nmonte={NMONTE}, obliquity={OBLIQUITY} deg)',
             fontsize=14, fontweight='bold')
fig.savefig(os.path.join(OUT_DIR, 'fig3_irradiance_4phases.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  -> fig3_irradiance_4phases.png saved")

# ═══════════════════════════════════════════
# FIGURE 4: Flux Time Series for Selected Faces
# ═══════════════════════════════════════════
print("\n[Figure 4] 繪製選定面元的照度時序圖...")
fig, ax = plt.subplots(figsize=(12, 5))

# Pick faces at different latitudes
face_lats = np.arcsin(face_centroids[:, 2] / np.linalg.norm(face_centroids, axis=1))
target_lats = [-60, -30, 0, 30, 60]
selected_faces = []
for lat in target_lats:
    lat_rad = np.radians(lat)
    idx = np.argmin(np.abs(face_lats - lat_rad))
    selected_faces.append(idx)

phases = np.linspace(0, 360, TT, endpoint=False)
colors = ['blue', 'cyan', 'green', 'orange', 'red']
for i, (nf, lat) in enumerate(zip(selected_faces, target_lats)):
    ax.plot(phases, Flux[nf, :], '-o', color=colors[i], markersize=3,
            label=f'Face {nf} (lat ~ {lat} deg)')

ax.set_xlabel('Orbital Phase (degrees)', fontsize=12)
ax.set_ylabel('Normalized Irradiance F', fontsize=12)
ax.set_title('Irradiance Time Series at Different Latitudes', fontsize=13)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 360)
fig.savefig(os.path.join(OUT_DIR, 'fig4_flux_timeseries.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  -> fig4_flux_timeseries.png saved")

# ═══════════════════════════════════════════
# FIGURE 5: Shadow Statistics
# ═══════════════════════════════════════════
print("\n[Figure 5] 繪製陰影統計分布...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Mean flux per face
mean_flux = np.mean(Flux, axis=1)
ax1.hist(mean_flux, bins=40, color='steelblue', edgecolor='white', alpha=0.8)
ax1.set_xlabel('Mean Irradiance')
ax1.set_ylabel('Number of Faces')
ax1.set_title(f'Distribution of Mean Irradiance\n(avg={mean_flux.mean():.3f})')
ax1.axvline(mean_flux.mean(), color='red', linestyle='--', label=f'mean={mean_flux.mean():.3f}')
ax1.legend()

# Fraction of time in shadow per face
shadow_frac = np.mean(Flux == 0, axis=1)
ax2.hist(shadow_frac, bins=30, color='dimgray', edgecolor='white', alpha=0.8)
ax2.set_xlabel('Fraction of Orbit in Shadow')
ax2.set_ylabel('Number of Faces')
ax2.set_title('Shadow Duration Distribution')
ax2.axvline(shadow_frac.mean(), color='red', linestyle='--', label=f'mean={shadow_frac.mean():.2f}')
ax2.legend()

fig.suptitle('Shadow & Irradiance Statistics', fontsize=14, fontweight='bold')
fig.savefig(os.path.join(OUT_DIR, 'fig5_shadow_statistics.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  -> fig5_shadow_statistics.png saved")

# ═══════════════════════════════════════════
# STEP 7: PDE — 表面熱傳導方程
# ═══════════════════════════════════════════
print("\n" + "=" * 60)
print("  STEP 7: PDE 生成 — 小行星表面熱傳導模型")
print("=" * 60)

print("""
=== 控制方程：一維徑向熱傳導 PDE ===

在每個表面面元 (facet) 下方，溫度 T(z,t) 滿足：

    rho * Cp * dT/dt = k * d²T/dz²          ... (1)

其中：
  z   = 徑向深度 (從表面向下為正)
  t   = 時間
  rho = 密度 (kg/m³)
  Cp  = 比熱容 (J/(kg·K))
  k   = 熱傳導率 (W/(m·K))

=== 表面邊界條件 (z = 0) ===

能量守恆：入射 = 輻射 + 傳導

    (1-A) * S * F(nf, t) = eps * sigma * T_s^4 - k * dT/dz|_{z=0}    ... (2)

其中：
  A     = 反照率 (Bond albedo ~ 0.1)
  S     = 太陽常數 (W/m²), 隨日心距 r_h 變化: S = 1361 / r_h²
  F(nf,t) = Monte Carlo 計算的照度因子 (本模擬輸出)
  eps   = 紅外發射率 (~0.9)
  sigma = Stefan-Boltzmann 常數 = 5.67e-8 W/(m²·K⁴)
  T_s   = 表面溫度 T(0, t)

=== 底部邊界條件 (z → ∞) ===

    dT/dz|_{z=∞} = 0                         ... (3)

(恆溫邊界，或零熱流)

=== 無量綱化 ===

定義熱慣量 Gamma = sqrt(k * rho * Cp)
定義穿透深度 l_s = sqrt(k * P / (pi * rho * Cp))

令 Z = z / l_s,  tau = t / P (P = 自轉週期)

PDE 變為：

    dT/dtau = (1/pi) * d²T/dZ²               ... (4)

表面邊界 (Z=0)：

    Gamma * sqrt(pi/P) * dT/dZ|_0 = (1-A)*S*F - eps*sigma*T^4   ... (5)

=== 數值離散 (Crank-Nicolson 格式) ===

空間: Z_j = j * dZ,  j = 0,1,...,N_z
時間: tau_n = n * dtau

    T_j^{n+1} - T_j^n     1    T_{j+1}^{n+1} - 2T_j^{n+1} + T_{j-1}^{n+1}
    ─────────────────── = ──── * ──────────────────────────────────────────────
          dtau             2pi                      dZ²

                                T_{j+1}^n - 2T_j^n + T_{j-1}^n
                         + ──── * ────────────────────────────────
                           2pi                   dZ²
""")

# ═══════════════════════════════════════════
# STEP 8: Solve PDE for representative faces
# ═══════════════════════════════════════════
print("[Step 8] 求解 PDE — 選定面元的溫度演化...")

# Physical constants
SIGMA = 5.67e-8       # Stefan-Boltzmann
ALBEDO = 0.1
EMISSIVITY = 0.9
S_SOLAR = 1361.0      # W/m² at 1 AU
RHO = 1500.0          # kg/m³ (typical asteroid)
CP = 600.0            # J/(kg·K)
K_COND = 0.01         # W/(m·K) (low thermal inertia)
P_ROT = 6 * 3600      # 6 hour rotation period (seconds)

Gamma = np.sqrt(K_COND * RHO * CP)
l_s = np.sqrt(K_COND * P_ROT / (np.pi * RHO * CP))

print(f"  熱慣量 Gamma = {Gamma:.2f} J/(m²·K·s^0.5)")
print(f"  穿透深度 l_s = {l_s*100:.2f} cm")

# 1D heat equation solver (Crank-Nicolson)
Nz = 30
Nt_per_orbit = TT
N_orbits = 5  # run multiple orbits to converge
dZ = 0.15
dtau = 1.0 / Nt_per_orbit

def solve_thermal_1d(flux_series, Nz=30, dZ=0.15, N_orbits=5):
    """Solve 1D heat equation for one face"""
    Nt = len(flux_series)
    dtau = 1.0 / Nt

    r = dtau / (2 * np.pi * dZ**2)

    T = np.ones(Nz) * 200.0  # initial guess 200K
    T_surface_history = []

    for orbit in range(N_orbits):
        for n in range(Nt):
            F_now = flux_series[n]
            T_new = T.copy()

            # Interior: Crank-Nicolson
            for j in range(1, Nz-1):
                T_new[j] = T[j] + r * (T[j+1] - 2*T[j] + T[j-1])

            # Bottom BC: dT/dz = 0
            T_new[Nz-1] = T_new[Nz-2]

            # Surface BC: energy balance (Newton iteration with damping)
            T_s = T_new[0]
            for _ in range(15):
                Q_in = (1 - ALBEDO) * S_SOLAR * F_now
                Q_rad = EMISSIVITY * SIGMA * T_s**4
                Q_cond = K_COND * (T_new[1] - T_s) / (l_s * dZ)
                residual = Q_in - Q_rad - Q_cond
                dres = -4 * EMISSIVITY * SIGMA * T_s**3 - K_COND / (l_s * dZ)
                if abs(dres) > 1e-20:
                    delta_T = -residual / dres
                    # Damping to prevent overshoot
                    delta_T = np.clip(delta_T, -50, 50)
                    T_s = T_s + delta_T
                T_s = np.clip(T_s, 40.0, 500.0)

            T_new[0] = T_s
            T = T_new

            if orbit == N_orbits - 1:
                T_surface_history.append(T_s)

    return np.array(T_surface_history), T

# Solve for selected faces
print("  正在求解 5 個緯度帶的溫度...")
T_results = {}
T_depth_profiles = {}
for i, (nf, lat) in enumerate(zip(selected_faces, target_lats)):
    T_hist, T_depth = solve_thermal_1d(Flux[nf, :])
    T_results[lat] = T_hist
    T_depth_profiles[lat] = T_depth
    T_max = T_hist.max()
    T_min = T_hist.min()
    print(f"  lat={lat:+3d} deg: T_min={T_min:.1f}K, T_max={T_max:.1f}K, Delta_T={T_max-T_min:.1f}K")

# ═══════════════════════════════════════════
# FIGURE 6: Surface Temperature Evolution
# ═══════════════════════════════════════════
print("\n[Figure 6] 繪製表面溫度隨軌道相位變化...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

phases = np.linspace(0, 360, TT, endpoint=False)
for i, lat in enumerate(target_lats):
    ax1.plot(phases, T_results[lat], '-', color=colors[i], linewidth=2,
             label=f'lat = {lat:+d} deg')

ax1.set_xlabel('Orbital Phase (degrees)', fontsize=12)
ax1.set_ylabel('Surface Temperature (K)', fontsize=12)
ax1.set_title('Surface Temperature vs Orbital Phase', fontsize=13)
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 360)

# Depth profiles at last timestep
depths = np.arange(Nz) * dZ * l_s * 100  # in cm
for i, lat in enumerate(target_lats):
    ax2.plot(T_depth_profiles[lat], depths, '-o', color=colors[i],
             markersize=3, label=f'lat = {lat:+d} deg')

ax2.set_xlabel('Temperature (K)', fontsize=12)
ax2.set_ylabel('Depth (cm)', fontsize=12)
ax2.set_title('Temperature Depth Profile (end of simulation)', fontsize=13)
ax2.invert_yaxis()
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

fig.suptitle('PDE Solution: 1D Heat Conduction', fontsize=14, fontweight='bold')
fig.savefig(os.path.join(OUT_DIR, 'fig6_temperature_pde.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  -> fig6_temperature_pde.png saved")

# ═══════════════════════════════════════════
# FIGURE 7: PDE Equations Summary
# ═══════════════════════════════════════════
print("\n[Figure 7] 繪製 PDE 公式彙整圖...")
fig, ax = plt.subplots(figsize=(14, 10))
ax.axis('off')

equations = [
    (0.95, r"$\bf{Monte\ Carlo\ Raytracing\ +\ Heat\ Conduction\ PDE\ Summary}$", 16),
    (0.88, r"$\bf{1.\ Governing\ PDE\ (1D\ Radial\ Heat\ Conduction):}$", 13),
    (0.83, r"$\rho\, C_p \,\frac{\partial T}{\partial t} = k \,\frac{\partial^2 T}{\partial z^2}$", 18),
    (0.76, r"$\bf{2.\ Surface\ Boundary\ Condition\ (z=0):}$", 13),
    (0.71, r"$(1-A)\, S_\odot\, F_{MC}(t) = \varepsilon\,\sigma\, T_s^4 \;-\; k\,\left.\frac{\partial T}{\partial z}\right|_{z=0}$", 16),
    (0.64, r"$\bf{3.\ Monte\ Carlo\ Irradiance\ Factor:}$", 13),
    (0.59, r"$F_{MC}(t) = \cos\theta_{solar} \;\times\; \left(1 - \frac{N_{occluded}}{N_{monte}}\right)$", 16),
    (0.52, r"$\bf{4.\ Shadow\ Detection\ Criterion:}$", 13),
    (0.47, r"$\hat{d}_{ij} \cdot \hat{s}_{solar} \geq \cos\alpha_{threshold}$" +
           r"$\quad (\alpha_{th} = $" + f"{np.degrees(np.arccos(initial_sc)):.1f}" + r"$^\circ)$", 14),
    (0.40, r"$\bf{5.\ Crater\ Profile\ (Power{-}Law\ Distribution):}$", 13),
    (0.35, r"$N(>D) \propto D^{-2}, \quad \delta_r = -d\,(1-x^2)\ (bowl),\quad d = 0.2\,D$", 14),
    (0.28, r"$\bf{6.\ Thermal\ Parameters:}$", 13),
    (0.23, r"$\Gamma = \sqrt{k\rho C_p}$" + f" = {Gamma:.2f}" +
           r"$\ \mathrm{J/(m^2 K\, s^{1/2})}$" +
           r"$,\quad l_s = \sqrt{\frac{kP}{\pi\rho C_p}}$" + f" = {l_s*100:.2f} cm", 13),
    (0.16, r"$\bf{7.\ Dimensionless\ PDE:}$", 13),
    (0.11, r"$\frac{\partial T}{\partial \tau} = \frac{1}{\pi}\,\frac{\partial^2 T}{\partial Z^2}$" +
           r"$,\quad Z = z/l_s,\quad \tau = t/P$", 16),
]

for y, text, size in equations:
    ax.text(0.05, y, text, transform=ax.transAxes, fontsize=size,
            verticalalignment='top', family='serif')

fig.savefig(os.path.join(OUT_DIR, 'fig7_pde_equations.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  -> fig7_pde_equations.png saved")

# ═══════════════════════════════════════════
# Save flux data
# ═══════════════════════════════════════════
flux_path = os.path.join(OUT_DIR, 'Flux_simulation.dat')
with open(flux_path, 'w') as f:
    f.write(f'{nface}\n')
    f.write(f'{TT}\n')
    for nf in range(nface):
        for t in range(TT):
            f.write(f'{Flux[nf, t]:.7f}\n')
print(f"\n  Flux data saved to Flux_simulation.dat")

print("\n" + "=" * 60)
print("  模擬完成！產出檔案：")
print("  1. fig1_asteroid_mesh.png      — 小行星 3D 網格")
print("  2. fig2_crater_profile.png     — 撞擊坑剖面 & 尺寸分布")
print("  3. fig3_irradiance_4phases.png — 四個軌道相位的照度圖")
print("  4. fig4_flux_timeseries.png    — 照度時序圖")
print("  5. fig5_shadow_statistics.png  — 陰影統計分布")
print("  6. fig6_temperature_pde.png    — PDE 溫度解")
print("  7. fig7_pde_equations.png      — PDE 公式彙整")
print("  8. Flux_simulation.dat         — 照度數據")
print("=" * 60)
