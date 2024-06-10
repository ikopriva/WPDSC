function [swpc, bc] = swpdec2a(x, n, varargin)
% SWPDEC2 Discrete stationary wavelet packet tree decomposition 2-D
%   SWPDEC2 performs a multilevel separable 2-D stationary wavelet packet
%   decomposition using a specific orthogonal wavelet ('wname' 
%   see WFILTERS for more information).
%
%   [SWPC,L] = SWPDEC2(X,N,'wname') computes full wavelet packet tree 
%   decomposition of image X at level N, using 'wname'.
%   N must be a strictly positive integer (see WMAXLEV for more
%   information).
%
%   The matrices of decomposition coefficients are stored in a 3D SWPC array.
%   Example: if n = 3, SWPC(:,:,1:4) correspond to level 1 decomposition, 
%   next 16 planes SWPC(:,:,5:20) correspond to level 2 and the last 64 planes
%   SWPC(:,:,21:84) correspond to level 3 decomposition. 
%   A, H, V, D coefficients are periodically repeating in the array planes.
%   Example: A = SWPC(:,:,1), H = SWPC(:,:,2), V = SWPC(:,:,3), D = SWPC(:,:,4); 
%   AA = SWPC(:,:,5), AH = SWPC(:,:,5), ... , DDV = SWPC(:,:,83), DDD = SWPC(:,:,84).
%
%   L is a book keeping vector, used for partial reconstruction. Non zero
%   elements means that corresponding SWPC planes are used for reconstruction.
%   Initially it is set to 0 except for the elements that represent the
%   first decomposition level.
%
%   See also SWT2, WPDEC2, SWPREC2.

%   D. Sersic 22-May-06.

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

l = 4.^[1:n]; % number of planes at each decomposition level
ind = cumsum(l)+1; ind = [1 max(1, ind(1:end-1))]; % first indices at each level
swpc = zeros(size(x,1), size(x,2), sum(l)); % space allocation

% level 1 decomposition
[swpc(:,:,1), swpc(:,:,2), swpc(:,:,3), swpc(:,:,4)] = swt2a(x, 1, lo_d, hi_d);
swpc(:,:,1:4)=swpc(:,:,1:4)/2; % Energy normalization

% levels 2 to n
for i = 2:n
    lo_d = dyadup(lo_d,0,1); hi_d = dyadup(hi_d,0,1);

    for j = 0:l(i-1)-1
        [swpc(:,:,ind(i)+4*j),swpc(:,:,ind(i)+4*j+1),swpc(:,:,ind(i)+4*j+2),swpc(:,:,ind(i)+4*j+3)] = ...
            swt2a(swpc(:,:,ind(i-1)+j),1,lo_d,hi_d);
        swpc(:,:,ind(i)+4*j:ind(i)+4*j+3)  =swpc(:,:,ind(i)+4*j:ind(i)+4*j+3)/2; % Energy normalization
    end
end

bc  = zeros(sum(l), 1); % book keeping vector
bc(1:4) = [1; 1; 1; 1];