function [swpc, bc] = swpdec(x, n, varargin)
% SWPDEC Discrete stationary wavelet packet tree decomposition 1-D
%   SWPDEC performs a multilevel 1-D stationary wavelet packet
%   decomposition using a specific orthogonal wavelet ('wname' 
%   see WFILTERS for more information).
%
%   [SWPC,L] = SWPDEC(X,N,'wname') computes full wavelet packet tree 
%   decomposition of signal X at level N, using 'wname'.
%   N must be a strictly positive integer (see WMAXLEV for more
%   information).
%
%   The vectors of decomposition coefficients are stored row-wise in SWPC matrix. 
%   Example: if n = 3, the first 2 rows correspond to level 1 decomposition, 
%   next 4 rows to level 2 and the last 8 rows correspond to
%   level 3 decomposition. Low and high pass filtered coefficients are alternating.
%   Rows SWPC(1,:) and SWPC(2,:) represent low and high pass level 1 coefficients.
%   Next 4 rows SWPC(3:6,:) represent level 2 decomposition, etc.
%
%   L is a book keeping vector, used for partial reconstruction. Non zero
%   elements means that corresponding SWPC rows are used for reconstruction.
%   Initially it is set to 0 except for the elements that represent the
%   first decomposition level.
%
%   See also SWT, WPDEC, SWPREC.

%   D. Sersic 21-Apr-06.
%   Rev. 22-May-06.

% Check arguments.
if nargin < 3
  error('Not enough input arguments.');
elseif nargin > 4
  error('Too many input arguments.');
end

% Compute decomposition filters.
if nargin==3
    [lo_d, hi_d] = wfilters(varargin{1},'d');
else
    lo_d = varargin{1}; hi_d = varargin{2};
end

l = 2.^[1:n]; % number of rows at each decomposition level
ind = l-1; % first element indices at each decomposition level
swpc = zeros(sum(l), length(x)); % space allocation

% level 1 decomposition
[swpc(1,:),swpc(2,:)] = swt(x, 1, lo_d, hi_d);

% levels 2 to n
for i = 2:n
    lo_d = dyadup(lo_d,0,1); hi_d = dyadup(hi_d,0,1);

    for j = 0:l(i-1)-1
        [swpc(ind(i)+2*j,:),swpc(ind(i)+2*j+1,:)] = swt(swpc(ind(i-1)+j,:), 1, lo_d, hi_d);
    end
end

bc = zeros(sum(l), 1); % book keeping vector
bc(1:2) = [1; 1];