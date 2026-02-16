function V = apply_triaxial(V, abc)
% APPLY_TRIAXIAL  Scale vertices by triaxial ellipsoid semi-axes.
%   V = apply_triaxial(V, [a b c])
%   Typical asteroid: a:b:c = 1.0:0.8:0.6

  V(:,1) = V(:,1) * abc(1);
  V(:,2) = V(:,2) * abc(2);
  V(:,3) = V(:,3) * abc(3);
end
