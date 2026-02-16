function V = add_accretion(V, n_lumps, amp_min, amp_max, sigma_min, sigma_max)
% ADD_ACCRETION  Add Gaussian bell-shaped accretion lumps.
%   V = add_accretion(V, n_lumps, amp_min, amp_max, sigma_min, sigma_max)
%   n_lumps   - number of lumps (default 10)
%   amp_min/max  - height range as fraction of mean radius (default 0.02..0.08)
%   sigma_min/max - angular width range in radians (default 0.2..0.6)

  if nargin < 2 || isempty(n_lumps),   n_lumps   = 10;   end
  if nargin < 3 || isempty(amp_min),   amp_min   = 0.02; end
  if nargin < 4 || isempty(amp_max),   amp_max   = 0.08; end
  if nargin < 5 || isempty(sigma_min), sigma_min = 0.2;  end
  if nargin < 6 || isempty(sigma_max), sigma_max = 0.6;  end

  R = mean(sqrt(sum(V.^2, 2)));

  % random centres on unit sphere
  centres = randn(n_lumps, 3);
  centres = centres ./ sqrt(sum(centres.^2, 2));

  % random amplitudes and widths
  amps   = amp_min + (amp_max - amp_min) * rand(n_lumps, 1);
  sigmas = sigma_min + (sigma_max - sigma_min) * rand(n_lumps, 1);

  for li = 1:n_lumps
    ctr   = centres(li, :);
    amp   = amps(li) * R;
    sigma = sigmas(li);

    for vi = 1:size(V, 1)
      p = V(vi, :);
      r_v = sqrt(sum(p.^2));
      p_hat = p / r_v;

      cos_ang = dot(p_hat, ctr);
      ang = acos(max(min(cos_ang, 1), -1));

      % Gaussian bell: delta = amp * exp(-ang^2 / (2*sigma^2))
      delta = amp * exp(-ang^2 / (2 * sigma^2));
      if delta > 1e-6 * R
        V(vi, :) = p_hat * (r_v + delta);
      end
    end
  end
end
