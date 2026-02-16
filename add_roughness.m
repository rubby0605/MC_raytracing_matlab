function V = add_roughness(V, F, amplitude, smooth_iter)
% ADD_ROUGHNESS  Add spatially-correlated surface roughness.
%   V = add_roughness(V, F, amplitude, smooth_iter)
%   amplitude   - noise amplitude as fraction of mean radius (default 0.005)
%   smooth_iter - Laplacian smoothing iterations (default 3)

  if nargin < 3 || isempty(amplitude),   amplitude   = 0.005; end
  if nargin < 4 || isempty(smooth_iter), smooth_iter = 3;     end

  nV = size(V, 1);
  R  = mean(sqrt(sum(V.^2, 2)));

  % random radial noise
  noise = amplitude * R * randn(nV, 1);

  % build adjacency for Laplacian smoothing
  % use sparse matrix: neighbours(i,j) = 1 if edge exists
  edges = [F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,1)];
  edges = [edges; edges(:,2) edges(:,1)];  % symmetric
  A = sparse(edges(:,1), edges(:,2), ones(size(edges,1),1), nV, nV);
  A = double(A > 0);  % binary adjacency

  % degree vector
  deg = full(sum(A, 2));

  % Laplacian smoothing of noise field
  for iter = 1:smooth_iter
    noise_new = zeros(nV, 1);
    for vi = 1:nV
      nbrs = find(A(vi, :));
      if ~isempty(nbrs)
        noise_new(vi) = mean(noise(nbrs));
      else
        noise_new(vi) = noise(vi);
      end
    end
    noise = noise_new;
  end

  % apply radial displacement
  for vi = 1:nV
    r_v = sqrt(sum(V(vi,:).^2));
    V(vi,:) = V(vi,:) / r_v * (r_v + noise(vi));
  end
end
