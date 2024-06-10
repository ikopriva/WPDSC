function [err_lrr23,Z,E] = nlrr_dualS23(X,lambda,display)
% This routine solves the following nuclear-norm optimization problem
% by using inexact Augmented Lagrange Multiplier, which has been also presented
% in the paper entitled "Robust Subspace Segmentation by Low-Rank Representation".
% min |Z|^{2/3}_{2/3}+lambda*|E|_2,1
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
% % rho = 1.1;                   % turnable
% % max_mu = 1e10;
% % mu = 1e-6;                   % turnable
% % xtx = X'*X;
% % inv_x = inv(xtx+eye(n));
% % %% Initializing optimization variables
% % for ii = 1:2
% %     if ii == 1
% %         M{ii} = rand(size(X, ii), 100);              % fix
% %         %M{ii} = randn(size(D, ii), ranks);
% %         [U{ii}, aa1] = qr(X'*M{ii}, 0);
% %     else
% %         M{ii} = rand(size(X, 3-ii), 100);            % fix
% %         %M{ii} = randn(size(D, 3-ii), ranks);
% %         [U{ii}, aa2] = qr(X'*M{ii}, 0);
% %     end
% % end
% % 
% % Z  = U{1}*U{2}';
% % a1 = norm(X, 2);
% % a2 = norm(X, Inf)/lambda;
% % Y1 = X/max(a1,a2);
% % 
% % %Z = zeros(n,n);
% % E = sparse(d,n);
% % 
% % %Y1 = zeros(d,n);
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
    ro = 2/mu;         % very key              
    if temp==0
    [U, S, V] = svd(temp+eps,'econ');
    else
    [U, S, V] = svd(temp,'econ');
    end
    a = diag(S);
    f = acosh(27*a.^2/16*ro^(-3/2));
    ei = 2/3^(1/2)*(ro)^(1/4)*(cosh(f/3)).^(1/2);
    b = ((ei+(2*a./ei-ei.^2).^(1/2))/2).^3;
    St = diag((a>2/3*(3*ro^3)^(1/4)).*b);
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
    err_lrr23(iter) = stopC;
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