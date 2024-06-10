function [err_lrr12,Z,E] = nlrr_dualS12(X,lambda,display)
% This routine solves the following nuclear-norm optimization problem
% by using inexact Augmented Lagrange Multiplier, which has been also presented
% in the paper entitled "Robust Subspace Segmentation by Low-Rank Representation".
% min |Z|^0.5_0.5+lambda*|E|_2,1
% s.t., X = XZ+E
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
% % if nargin<2
% %     lambda = 1;
% % end
% % if nargin<3
% %     display = true;
% % end
% % tol = 1e-8;
% % maxIter = 1e6;
% % [d, n] = size(X);
% % rho = 1.1;             % turnable
% % max_mu = 1e10;
% % mu = 1e-6;             % turnable
% % xtx = X'*X;
% % inv_x = inv(xtx+eye(n));
% % %% Initializing optimization variables
% % Z = zeros(n,n);
% % E = sparse(d,n);
% % 
% % Y1 = zeros(d,n);
% % Y2 = zeros(n,n);

if nargin<2
    lambda = 1;
end
if nargin<3
    display = true;
end
tol = 1e-8;
maxIter = 1e6;
[d, n] = size(X);
rho = 1.1;            % turnable
max_mu = 1e10;
mu = 1e-6;            % turnable
xtx = X'*X;
inv_x = inv(xtx+eye(n));
%% Initializing optimization variables
Z = zeros(n,n);
E = sparse(d,n);

Y1 = X;%zeros(d,n);
Y2 = zeros(n,n);
%% Start main loop
iter = 0;
% disp(['initial,rank=' num2str(rank(Z))]);
while iter<maxIter
    iter = iter + 1;
    
    % update J
    temp = Z + Y2/mu;
    ro = 1/mu;
    if temp==0
    [U, S, V] = svd(temp+eps,'econ');
    else
    [U, S, V] = svd(temp,'econ');
    end
    a = diag(S);
    b = 2/3*a.*(1+cos(2*pi/3-2/3*acos(ro/4*(a/3).^(-3/2))));
    St = diag((a>(54^(1/3)/4)*(2*ro)^(2/3)).*b);
    J = U(:,1:size(St,1))*(St*V(:,1:size(St,1))');
    
    % update Z
    Z = inv_x*(xtx-X'*E+J+(X'*Y1-Y2)/mu);
    
    % update E
    xmaz = X-X*Z;
    temp = xmaz+Y1/mu;
    E = solve_l1l2(temp,lambda/mu);
    %
    leq1 = xmaz-E;
    leq2 = Z-J;
    %stopC = max(max(max(abs(leq1))),max(max(abs(leq2))));
    stopC = max(max(abs(leq1)));
    err_lrr12(iter) = stopC;
    if ~display && (iter==1 || mod(iter,300)==0 || stopC<tol)
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ...
            ',rank=' num2str(rank(Z,1e-3*norm(Z,2))) ',stopALM=' num2str(stopC,'%2.3e')]);
    end
    if stopC<tol
        break;
    else
        Y1 = Y1 + mu*leq1;
        Y2 = Y2 + mu*leq2;
        mu = min(max_mu,mu*rho);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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