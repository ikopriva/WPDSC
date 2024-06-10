function [Z,E] = nlrr_S12(X,lambda,rho,display,DEBUG)
% This matlab code implements linearized ADM method for LRR problem
% min |Z|^p_Sp+lambda*|E|_2,1
% s.t., X = XZ+E
% inputs:
%        X -- D*N data matrix
% outputs:
%        Z -- N*N representation matrix
%        E -- D*N sparse error matrix
%        relChgs --- relative changes
%        recErrs --- reconstruction errors
[d, n] = size(X);
if nargin < 5
    DEBUG = true;
end
if nargin < 4
    display = true;
end
if nargin < 3
    rho = 1.1;
end
if nargin < 2
    lambda = 1/sqrt(log(n));
end

norm2X = norm(X,2);
tol = 1e-8;
maxIter = 1000;
max_mu = 1e10;

% choice of both mu and eta
mu = 1e-2;
eta = norm2X*norm2X*1.05; %eta needs to be larger than ||X||_2^2, but need not be too large.

%% Initializing optimization variables
E = sparse(d,n);
Y = zeros(d,n);
Z = zeros(n,n);
XZ = zeros(d,n); % XZ = X*Z;

%% Start main loop
iter = 0;

while iter<maxIter
    iter = iter + 1;
    
    % update Z based on S_1/2 norm
    ro = 1/(mu*eta);
    M = Z + X'*(X - XZ - E + Y/mu)/eta;
    [U, S, V] = svd(M,'econ');
    a = diag(S);
    b = 2/3*a.*(1+cos(2*pi/3-2/3*acos(ro/4*(a/3).^(-3/2))));
    St = diag((a>(54^(1/3)/4)*(2*ro)^(2/3)).*b);
    Z = U(:,1:size(St,1))*(St*V(:,1:size(St,1))');
    
    % update E
    XZ = X*Z;                    % introducing XZ to avoid computing X*Z multiple times, which has O(n^3) complexity.
    E = solve_l1l2(X - XZ + Y/mu,lambda/mu);
    
    % stop criteria
    dY = X - XZ - E;
    sc(iter) = max(max(abs(dY)));
    convergenced = sc(iter) < tol;
    if convergenced
        break;
    else
        Y = Y + mu*dY;
        mu = min(max_mu, mu*rho);
    end
    
    if DEBUG
        if ~display && (iter==1 || mod(iter,100)==0 || convergenced)
            disp([' iter ' num2str(iter) ', mu=' num2str(mu) ...
                ', rank(Z)=' num2str(rank(Z,1e-3*norm(Z,2))) ', stopC='  num2str(sc(iter))]);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end
function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end