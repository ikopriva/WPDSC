function x = swprec(swpc, bc, varargin)
% SWPREC Discrete stationary wavelet packet tree reconstruction 1-D
%   X = SWPREC(SWPC,L,'wname') or X = SWPREC(SWPC,L,Lo_R,Hi_R) 
%   reconstructs the signal X based on the multilevel stationary wavelet   
%   packet decomposition structure SWPC (see SWPDEC).
%
%   L is a book keeping vector, used for partial reconstruction. 
%   Non zero elements mean that corresponding SWPC rows are used for
%   the reconstruction, otherwise they are not used.
%
%   See also SWPDEC, WPREC.

%   D. Sersic 23-May-06.

    % Check arguments.
    if nargin < 3
      error('Not enough input arguments.');
    elseif nargin > 4
      error('Too many input arguments.');
    end

    % Compute reconstruction filters.
    if nargin==3
        [lo_r, hi_r] = wfilters(varargin{1},'r');
    else
        lo_r = varargin{1}; hi_r = varargin{2};
    end

    % Check the book keeping vector: if zero go recursively to the next level
    for i=1:2
        if 0 == bc(i)
            [bcr, swpcr] = reduce_order(bc, swpc, i-1);
            if ~isempty(bcr) % next level
                lo_rr = dyadup(lo_r,1,1);  hi_rr = dyadup(hi_r,1,1);
                swpc(i,:) = swprec(swpcr, bcr, lo_rr, hi_rr);
            else    % end of the branch
                swpc(i,:) = 0;
            end
        end
    end
    % The reconstruction
    x = iswt(swpc(1,:),swpc(2,:), lo_r, hi_r);
end

function [bcr, swpcr] = reduce_order(bc, swpc, offset)
    [r,c] = size(bc);

    n = floor(log2(r)); % number of decomposition levels
    l = 2.^[1:n]; % number of rows at each decomposition level
    ind = l-1; % first element indices at each decomposition level

    index = [];
    for i = 2:n
        index = [index, l(i-1)*offset + (ind(i):ind(i)+l(i-1)-1)];
    end

    bcr   =   bc(index, :);
    swpcr = swpc(index, :);
end
