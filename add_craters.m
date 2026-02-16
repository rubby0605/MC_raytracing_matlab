function V = add_craters(V, F, n_craters, d_min, d_max, depth_ratio)
% ADD_CRATERS  Add impact craters with power-law size distribution.
%   V = add_craters(V, F, n_craters, d_min, d_max, depth_ratio)
%   n_craters  - number of craters (default 50)
%   d_min      - minimum diameter as fraction of mean radius (default 0.05)
%   d_max      - maximum diameter as fraction of mean radius (default 0.4)
%   depth_ratio - depth/diameter ratio (default 0.2, simple crater)
%
%   Size distribution: N(>D) ~ D^(-2)  =>  CDF F(D) = 1-(D_min/D)^2
%   Inverse CDF: D = D_min / sqrt(1 - U),  U ~ Uniform(0, 1-d_min^2/d_max^2)

  if nargin < 3 || isempty(n_craters), n_craters = 50;   end
  if nargin < 4 || isempty(d_min),     d_min = 0.05;     end
  if nargin < 5 || isempty(d_max),     d_max = 0.4;      end
  if nargin < 6 || isempty(depth_ratio), depth_ratio = 0.2; end

  % mean body radius
  R = mean(sqrt(sum(V.^2, 2)));

  % generate crater diameters via inverse CDF of power-law
  u_max = 1 - (d_min / d_max)^2;
  u = rand(n_craters, 1) * u_max;
  diameters = d_min ./ sqrt(1 - u) * R;

  % sort large to small (large craters first)
  [diameters, idx] = sort(diameters, 'descend');

  % random crater centres on unit sphere
  centres = randn(n_craters, 3);
  centres = centres ./ sqrt(sum(centres.^2, 2));

  % apply each crater
  for ci = 1:n_craters
    D   = diameters(ci);
    rad = D / 2;                    % crater radius
    dep = depth_ratio * D;          % crater depth
    rim_h = dep * 0.15;             % rim height (~15% of depth)
    ctr = centres(ci, :);

    for vi = 1:size(V, 1)
      p = V(vi, :);
      r_v = sqrt(sum(p.^2));
      p_hat = p / r_v;

      % angular distance from crater centre
      cos_ang = dot(p_hat, ctr);
      ang = acos(max(min(cos_ang, 1), -1));
      dist = ang * R;  % arc-length distance

      if dist < rad * 1.5
        % normalised distance 0..1 within crater bowl
        x = dist / rad;
        if x <= 1.0
          % parabolic bowl: delta = -dep * (1 - x^2)
          delta = -dep * (1 - x * x);
        elseif x <= 1.5
          % rim: raised cosine taper
          t = (x - 1.0) / 0.5;  % 0..1
          delta = rim_h * 0.5 * (1 + cos(pi * t));
        else
          delta = 0;
        end
        % apply radial displacement
        V(vi, :) = p_hat * (r_v + delta);
      end
    end
  end
end
