function [Z,E] = lrr_S1(X,lambda,rho,display,DEBUG)
% This matlab code implements linearized ADM method for LRR problem
% min |Z|_*+lambda*|E|_2,1 s.t., X = XZ+E
% inputs:
%        X -- D*N data matrix
% outputs:
%        Z -- N*N representation matrix
%        E -- D*N sparse error matrix
if nargin < 5
    DEBUG = true;
end
if nargin < 4
    display = true;
end
if nargin < 3
    rho = 1.9;
end
if nargin < 2
    lambda = 0.1;
end

tol = 1e-8;                %threshold for the error in constraint
max_mu = 1e10;
maxIter = 1000;

% values of mu and eta
[d, n] = size(X);
norm2X = norm(X,2);
mu = 1e-6;
eta = norm2X*norm2X*1.02; %eta needs to be larger than ||X||_2^2, but need not be too large.

%% Initialization(very key)
E = sparse(d,n);
Y = zeros(d,n);
Z = zeros(n, n);
XZ = zeros(d, n);

%% Start main loop
iter = 0;
while iter<maxIter
    iter = iter + 1;
    
    % update Z
    M = Z + X'*(X - XZ - E + Y/mu)/eta;
    [U, S, V] = svd(M,'econ');
    S = diag(S);
    svp = length(find(S>1/(mu*eta)));
    if svp>=1
        S = S(1:svp)-1/(mu*eta);
    else
        svp = 1;
        S = 0;
    end
    A.U = U(:, 1:svp);
    A.s = S;
    A.V = V(:, 1:svp);
    Z = A.U*diag(A.s)*A.V';
    
    % update E
    XZ = X*Z;   % introducing XZ to avoid computing X*Z multiple times, which has O(n^3) complexity.
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
        if display && (iter==1 || mod(iter,50)==0 || convergenced)
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