function write_property_dat(number, a, b, c, obliquity, name, filename)
% WRITE_PROPERTY_DAT  Generate Property.dat for Shadow_maker compatibility.
%   write_property_dat(number, a, b, c, obliquity, name, filename)
%   Matches the read format at Shadow script lines 58-72:
%     Line 1: total_number (integer)
%     Line 2: header
%     Line 3+: Number a b c spare obliquity Name
%
%   number    - asteroid catalogue number (default 74971)
%   a, b, c   - semi-axes in km (default 1.0, 0.8, 0.6)
%   obliquity - obliquity in degrees (default 45)
%   name      - asteroid name string (default 'Asteroid')
%   filename  - output file (default 'Property.dat')

  if nargin < 1 || isempty(number),    number    = 74971;          end
  if nargin < 2 || isempty(a),         a         = 1.0;            end
  if nargin < 3 || isempty(b),         b         = 0.8;            end
  if nargin < 4 || isempty(c),         c         = 0.6;            end
  if nargin < 5 || isempty(obliquity), obliquity = 45;             end
  if nargin < 6 || isempty(name),      name      = 'Asteroid';     end
  if nargin < 7 || isempty(filename),  filename  = 'Property.dat'; end

  fid = fopen(filename, 'w');
  % Line 1: total number of entries
  fprintf(fid, '%d\n', 1);
  % Line 2: header
  fprintf(fid, 'Number a b c spare obliquity Name\n');
  % Line 3: data — spare field set to 0
  fprintf(fid, '%d %f %f %f %f %f %s\n', number, a, b, c, 0.0, obliquity, name);
  fclose(fid);
end
