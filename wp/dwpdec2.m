function [dwpc, bc] = dwpdec2(x, n, varargin)
% DWPDEC2 Discrete wavelet packet tree decomposition 2-D
%   DWPDEC2 performs a multilevel separable 2-D decimated wavelet packet
%   decomposition using a specific orthogonal wavelet ('wname' 
%   see WFILTERS for more information).
%
%   [DWPC,L] = DWPDEC2(X,N,'wname') computes full wavelet packet tree 
%   decomposition of image X at level N, using 'wname'.
%   N must be a strictly positive integer (see WMAXLEV for more
%   information).
%
%   The matrices of decomposition coefficients are stored in a 3D DWPC array.
%   Example: if n = 2, SWPC(:,:,1) is organized as [A H; V D] and corresponds 
%   to level 1 decomposition; SWPC(:,:,2) is organized as:
%     [AA HA AH HH 
%      VA DA VH DH
%      AV HV AD HD
%      VV DV VD DD]
%   and corresponds to level 2 decopmosition.
%
%   L is a book keeping vector, used for partial reconstruction. Non zero elements 
%   means that corresponding DWPC coefficients are used for the reconstruction.
%   Example: if n = 2, the correspondence between elements of L and
%   decomposition coefficients is:
%                  [5  6  9 10
%       [1 2        7  8 11 12
%        3 4]      13 14 17 18
%                  15 16 19 20]
%   Initially L is set to 0 except for the elements that represent the
%   first decomposition level (1,2,3,4).
%
%   See also SWPDEC2, DWT2.

%   D. Sersic 05-Sep-06.

% Check arguments.
if nargin < 3
  error('Not enough input arguments.');
elseif nargin > 4
  error('Too many input arguments.');
end

if n > 8
    error('More then 8 levels are not supported.');
end

% Preserve initial size.
s = size(x);
pow = 2^n;
if any(rem(s,pow))
    sOK = ceil(s/pow)*pow;
    oriStr = ['(' int2str(s(1))   ',' int2str(s(2)) ')'];
    sugStr = ['(' int2str(sOK(1)) ',' int2str(sOK(2)) ')'];
    msg = strvcat(...
        ['The level of decomposition ' int2str(n)],...
        ['and the size of the image ' oriStr],...
        'are not compatible.',...
        ['Suggested size: ' sugStr],...
        '(see Image Extension Tool)', ...
        ' ', ...
        ['2^Level has to divide the size of the image.'] ...
            );
    errargt(mfilename,msg,'msg');
    varargout = {[] };
    return
end

% Compute decomposition filters.
if nargin==3
    [lo_d, hi_d] = wfilters(varargin{1},'d');
else
    lo_d = varargin{1}; hi_d = varargin{2};
end

l = 4.^[1:n]; % number of images at each decomposition level
r = s(1); c = s(2);
    
dwpc = zeros(r, c, n); % space allocation

% level 1 decomposition
dwpc(:,:,1) = mdwt2(x, lo_d, hi_d);

% levels 2 to n
for i = 2:n
    r = r/2; c = c/2;
    for j = 0:l(i-1)-1
        [m, k] = j2km(j);
        dwpc(r*m+(1:r),c*k+(1:c),i) = mdwt2(dwpc(r*m+(1:r),c*k+(1:c),i-1),lo_d,hi_d);
    end
end

bc  = zeros(sum(l), 1); % book keeping vector
bc(1:4) = [1; 1; 1; 1];

end


function m = mdwt2(x, lo_d, hi_d)
% Result of decimated DWT2 is reduced to the quarter size images and
% concatenated in a single image
 [a, h, v, d] = dwt2(x, lo_d, hi_d);
 s = size(x)/2;
 a = wkeep2(a, s); h = wkeep2(h, s); 
 v = wkeep2(v, s); d = wkeep2(d, s);
 m = [a h; v d];
end
