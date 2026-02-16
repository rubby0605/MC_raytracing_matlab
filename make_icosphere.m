function [V, F] = make_icosphere(level)
% MAKE_ICOSPHERE  Generate an icosphere by recursive subdivision.
%   [V, F] = make_icosphere(level)
%   level = number of subdivision iterations (4 -> 2562 verts, 5 -> 10242)
%   V = Nx3 vertices on the unit sphere
%   F = Mx3 face indices (1-based)

  % --- Golden ratio icosahedron ---
  t = (1 + sqrt(5)) / 2;
  V = [ -1  t  0;  1  t  0; -1 -t  0;  1 -t  0;
         0 -1  t;  0  1  t;  0 -1 -t;  0  1 -t;
         t  0 -1;  t  0  1; -t  0 -1; -t  0  1 ];
  V = V ./ sqrt(sum(V.^2, 2));

  F = [  1 12  6;  1  6  2;  1  2  8;  1  8 11;  1 11 12;
         2  6 10;  6 12  5; 12 11  3; 11  8  7;  8  2  9;
         4 10  5;  4  5  3;  4  3  7;  4  7  9;  4  9 10;
         5 10  6;  3  5 12;  7  3 11;  9  7  8; 10  9  2 ];

  % --- Recursive subdivision (vectorized) ---
  for iter = 1:level
    nV = size(V, 1);
    nF = size(F, 1);

    % Extract all 3 edges per face: (v1,v2), (v2,v3), (v3,v1)
    e1 = [F(:,1) F(:,2)];
    e2 = [F(:,2) F(:,3)];
    e3 = [F(:,3) F(:,1)];
    edges = [e1; e2; e3];  % 3*nF x 2

    % Canonical ordering: lo < hi
    edges = [min(edges, [], 2), max(edges, [], 2)];

    % Find unique edges and map back
    [uniq_edges, ~, edge_idx] = unique(edges, 'rows');
    nE = size(uniq_edges, 1);

    % Create midpoint vertices (all at once)
    midpts = (V(uniq_edges(:,1), :) + V(uniq_edges(:,2), :)) / 2;
    midpts = midpts ./ sqrt(sum(midpts.^2, 2));

    % Append to vertex list; midpoint indices start at nV+1
    V = [V; midpts];

    % Map edge indices to new vertex indices
    mid_idx = nV + edge_idx;  % 3*nF x 1
    a = mid_idx(1:nF);        % midpoint of edge (v1,v2)
    b = mid_idx(nF+1:2*nF);   % midpoint of edge (v2,v3)
    c = mid_idx(2*nF+1:3*nF); % midpoint of edge (v3,v1)

    % Build 4 sub-triangles per original face
    v1 = F(:,1); v2 = F(:,2); v3 = F(:,3);
    F = [v1 a c; v2 b a; v3 c b; a b c];
  end
end
