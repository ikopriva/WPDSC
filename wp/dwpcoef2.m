function coeff = dwpcoef2(dwpc, j)
% DWPCOEF2 Get discrete wavelet packet tree 2-D coefficients
%   DWPCOEF2(DWPC, L) gets 2-D coefficients determined by linear index L.
%
%   Example: if n = 2, the correspondence between elements of L and
%   decomposition coefficients is:
%                  [5  6  9 10
%       [1 2        7  8 11 12
%        3 4]      13 14 17 18
%                  15 16 19 20]
%
%   See also SWPDEC2.

%   D. Sersic 07-Sep-06.

% Check arguments.
if nargin ~= 2
  error('Two input arguments needed.');
end

[r, c, n] = size(dwpc);
l = 4.^[1:n]; % number of images at each decomposition level

% Determine the decomposition level
L = cumsum(l); 
i = find(j<=L, 1);
if i-1, 
    j = j-L(i-1)-1; 
else
    j = j-1; 
end

r = r/2^i; c = c/2^i;  % dimension of the image
[m, k] = j2km(j); % image indices

% Extract the image (WP coefficients)
coeff = dwpc(r*m+(1:r),c*k+(1:c),i);
