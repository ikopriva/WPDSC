function x = swprec2(swpc, bc, varargin)
% SWPREC2 Discrete stationary wavelet packet tree reconstruction 2-D
%   X = SWPREC2(SWPC,L,'wname') or X = SWPREC2(SWPC,L,Lo_R,Hi_R) 
%   reconstructs the image X based on the multilevel stationary wavelet   
%   packet decomposition structure SWPC (see SWPDEC2).
%
%   L is a book keeping vector, used for partial reconstruction. 
%   Non zero elements mean that corresponding SWPC planes are used for
%   the reconstruction, otherwise they are not used.
%
%   See also SWPDEC2, WPREC2, ISWT2.

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
    for i=1:4
        if 0 == bc(i)
            [bcr, swpcr] = reduce_order(bc, swpc, i-1);
            if ~isempty(bcr) % next level
                lo_rr = dyadup(lo_r,1,1);  hi_rr = dyadup(hi_r,1,1);
                swpc(:,:,i) = swprec2(swpcr, bcr, lo_rr, hi_rr);
            else    % end of the branch
                swpc(:,:,i) = 0;
            end
        end
    end
    % The reconstruction
    x = myiswt2(swpc(:,:,1),swpc(:,:,2),swpc(:,:,3),swpc(:,:,4), lo_r, hi_r);
end

function [bcr, swpcr] = reduce_order(bc, swpc, offset)
    [r,c] = size(bc);

    n = floor(log2(r)/2); % number of decomposition levels
    l = 4.^[1:n]; % number of planes at each decomposition level
    ind = cumsum(l)+1; ind = [1 max(1, ind(1:end-1))]; % first indices at each level

    index = [];
    for i = 2:n
        index = [index, l(i-1)*offset + (ind(i):ind(i)+l(i-1)-1)];
    end

    bcr   =   bc(index,  :);
    swpcr = swpc(:,:,index);
end


function x = myiswt2(a, h, v, d, lo_r, hi_r)
    lf = length(lo_r);
    s  = size(a);

    a  = wextend('2D', 'per', a, [lf/2, lf/2]); h = wextend('2D', 'per', h, [lf/2, lf/2]); 
    v  = wextend('2D', 'per', v, [lf/2, lf/2]); d = wextend('2D', 'per', d, [lf/2, lf/2]); 

    t0 = conv2(a, lo_r') + conv2(h, hi_r');   d0 = conv2(v, lo_r') + conv2(d, hi_r');
    x  = conv2(t0, lo_r) + conv2(d0, hi_r);

    x = wkeep(x, s)/4;
end
