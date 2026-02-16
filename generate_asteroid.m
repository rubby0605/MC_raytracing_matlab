%% generate_asteroid.m — Procedural asteroid mesh generator
%  Outputs: spot, facet, npoint, nface, str  (compatible with Shadow_maker)
%  Also generates Property.dat
clear; close all;

%% ===== Parameters =====
seed          = 42;              % RNG seed for reproducibility
subdiv_level  = 4;               % icosphere subdivision (4->2562v, 5->10242v)
triaxial      = [1.0 0.8 0.6];  % semi-axes a:b:c

% Crater parameters
n_craters   = 50;     % number of impact craters
crater_dmin = 0.05;   % min diameter (fraction of mean radius)
crater_dmax = 0.4;    % max diameter (fraction of mean radius)
depth_ratio = 0.2;    % depth/diameter (simple crater)

% Accretion lump parameters
n_lumps   = 10;
lump_amp  = [0.02 0.08];  % [min max] height (fraction of R)
lump_sig  = [0.2  0.6];   % [min max] angular width (radians)

% Roughness parameters
roughness_amp   = 0.005;  % amplitude (fraction of R)
roughness_smooth = 3;     % Laplacian smoothing iterations

% Property.dat parameters
ast_number    = 74971;
ast_obliquity = 45;        % degrees
ast_name      = 'Asteroid';

%% ===== Pipeline =====
rng(seed);

% Step 1: Base icosphere
fprintf('Step 1: Generating icosphere (level %d)...\n', subdiv_level);
[V, F] = make_icosphere(subdiv_level);
fprintf('  -> %d vertices, %d faces\n', size(V,1), size(F,1));

% Step 2: Triaxial ellipsoid deformation
fprintf('Step 2: Applying triaxial scaling [%.1f %.1f %.1f]...\n', triaxial);
V = apply_triaxial(V, triaxial);

% Step 3: Impact craters
fprintf('Step 3: Adding %d craters...\n', n_craters);
V = add_craters(V, F, n_craters, crater_dmin, crater_dmax, depth_ratio);

% Step 4: Accretion lumps
fprintf('Step 4: Adding %d accretion lumps...\n', n_lumps);
V = add_accretion(V, n_lumps, lump_amp(1), lump_amp(2), lump_sig(1), lump_sig(2));

% Step 5: Surface roughness
fprintf('Step 5: Adding surface roughness...\n');
V = add_roughness(V, F, roughness_amp, roughness_smooth);

% Step 6: Fix face winding (ensure outward normals)
fprintf('Step 6: Checking face winding...\n');
n_flipped = 0;
for fi = 1:size(F, 1)
  a = V(F(fi,1), :);
  b = V(F(fi,2), :);
  c = V(F(fi,3), :);
  centroid = (a + b + c) / 3;
  normal = cross(c - a, b - a);
  if dot(normal, centroid) < 0
    F(fi, :) = [F(fi,1) F(fi,3) F(fi,2)];
    n_flipped = n_flipped + 1;
  end
end
fprintf('  -> Flipped %d faces\n', n_flipped);

% Step 7: Write Property.dat
fprintf('Step 7: Writing Property.dat...\n');
write_property_dat(ast_number, triaxial(1), triaxial(2), triaxial(3), ...
                   ast_obliquity, ast_name);

%% ===== Export to workspace variables =====
spot   = V;
facet  = F;
npoint = size(V, 1);
nface  = size(F, 1);
str    = ast_name;

fprintf('\nDone! Workspace variables:\n');
fprintf('  spot:   %d x 3\n', npoint);
fprintf('  facet:  %d x 3\n', nface);
fprintf('  npoint: %d\n', npoint);
fprintf('  nface:  %d\n', nface);
fprintf('  str:    %s\n', str);

%% ===== 3D Preview =====
figure('Name', 'Asteroid Preview', 'NumberTitle', 'off');
p = patch('Faces', facet, 'Vertices', spot, 'FaceColor', [0.6 0.55 0.5], ...
          'EdgeColor', 'none', 'FaceLighting', 'gouraud');
light('Position', [1 1 1]);
light('Position', [-1 -0.5 -0.5], 'Color', [0.3 0.3 0.3]);
axis equal; axis off;
material dull;
title(sprintf('Asteroid: %d verts, %d faces', npoint, nface));
view(30, 20);
fprintf('\nReady to run: run(''1_Shadow_maker(NObug0720).m'')\n');
