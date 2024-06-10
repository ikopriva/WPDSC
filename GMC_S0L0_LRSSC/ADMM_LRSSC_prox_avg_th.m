% Non-convex Low Rank Sparse Subbspace Clustering for estimation of the
% matrix of coefficients:
% min\|X-XC\|_F^2 + \tau\|S\|_1 + \lambda\sum(fi(sig(L),a)
%
% INPUTS:
%   X: PxN data matrix with n samples and p features
%   lambda: regularization constant related to sparsity induced regualrizer of singular values
%   tau: weight of L1 norm regularization for entries of S
%   affine: 1 - affine subspace model; 0 - independent subspace model
%   opts:  Structure value with following fields:
%          opts.lambda:    coefficients for low-rank constraint
%          opts.rank_est:    parametar for nonconvex rank estimator, if <1
%          convex problem
%          opts.mu1:  penalty parameter for auxiliary variable C1 in augmented Lagrangian
%          opts.mu2:  penalty parameter for auxiliary variable C2 in augmented Lagrangian
%          opts.mu3:  penalty parameter for affine subspaces constraint in augmented Lagrangian
%          opts.max_mu1:  maximum  penalty parameter for mu1 parameter
%          opts.max_mu2:  maximum  penalty parameter for mu2 parameter
%          opts.max_mu3:  maximum  penalty parameter for mu3 parameter (used only for affine subspaces)
%          opts.rho1: step size for adaptively changing mu1, if 1 fixed mu1 is used
%          opts.rho2: step size for adaptively changing mu2, if 1 fixed mu2 is used
%          opts.rho3: step size for adaptively changing mu3, if 1 fixed mu3
%          is used (used only for affine subspaces)
%          opts.error_bound: error bound for convergence
%          opts.iter_max:  maximal number of iterations
%          opts.affine: true for the affine subspaces (default: false)
%          opts.soft_thr:  true for soft thresholding (minimization of L1
%          norm, convex problem), false for hard threshold (minimization of L0 norm, non-convex problem)
%          (default: false)
%          opts.d:  subspace dimension used for threholding values of C
%                   matrix
%
% OUTPUTS:
%   C: NxN matrix of coefficients
%   RMSE: error
%   error: ||X-XC||/||X||
%
% Ivica Kopriva, March 2023.

function [C, error] = ADMM_LRSSC_prox_avg_th (X,opts)

if ~exist('opts', 'var')
    opts = [];
end

% default parameters
rho = 3;
max_mu1 = 1e6;
gamma = 0.6; % GMC parameter
err_thr = 1e-4;
iter_max = 100;

if isfield(opts, 'lambda');      lambda = opts.lambda;      end
if isfield(opts, 'alpha');      alpha = opts.alpha;      end % elra parameter
if isfield(opts, 'double_thr');    double_thr = opts.double_thr;    end
if isfield(opts, 'l0norm');      l0norm = opts.l0norm;      end
if isfield(opts, 'gmc');      gmc = opts.gmc;      end
if isfield(opts, 'iter_max');    iter_max = opts.iter_max;    end
if isfield(opts, 'err_thr');    err_thr = opts.err_thr;    end
if isfield(opts, 'gamma');      gamma = opts.gamma;      end
if isfield(opts, 'rho');      rho = opts.rho;      end
if isfield(opts,'d');         d = opts.d;   end

%% initialization

[M,N]=size(X);

%lambda1 = 1/(1+lambda);
%lambda2 = lambda/(1+lambda);

lambda1 = lambda;
lambda2 = 1-lambda;

% setting penalty parameters for the ADMM
mu1 = alpha;

J = zeros(N,N);  % auxiliary variable
%J = rand(N, N);
C = J;

% Lagrange multpliers
LAM_1 = zeros(N,N);

% Fixed precomputed term for J
tic;

XT = X'*X;

Jf = inv(XT + mu1*eye(N));
J = Jf*(XT + mu1*C - LAM_1);
J = normc(J);  % necessary if X is column normalized (YaleB dataset !!!!!!!)

not_converged = 1;
iter=1;

L_grad_vect = [];
err1_all =  [];
err2_all = [];

while not_converged
    
    J_prev = J;
    
    % Update of J
    J = Jf*(XT + mu1*C - LAM_1);
    J = normc(J);  % necessary if X is column normalized (YaleB dataset !!!!!!!)
    
    % Update of C1
    %[U,Sig,V]=svdsecon(J+LAM_1/mu1,k);
    [U Sig V] = svd(J+LAM_1/mu1,'econ');
    sig = diag(Sig)';
    thr = lambda1/mu1;
    
    if l0norm
        thr = sqrt(2*thr);
        if double_thr
            thr = thr*2;
        end
        sig_thr = sig.*((sign(abs(sig)-thr)+1)/2);
        
        [is inds] = sort(sig_thr,'descend');
        ind = inds(1:sum(sign(is)));
        sig = sig_thr(ind);
        V = V(:,ind);
        U = U(:,ind);
        Sig = diag(sig);
        C1 = U*Sig*V';
        
    elseif gmc
        tmp = arrayfun(@(y) firm_thresh(y, thr, thr/gamma), sig);
        C1 = U*diag(tmp)*V';
    end
    
        
    [cm indmax] = maxk(abs(C1),d,1);
    tmp = zeros(N,N);
    
    for n=1:N
        tmp(indmax(:,n),n)=C1(indmax(:,n),n);
    end
    C1 = tmp; clear tmp;
    
    % Update of C2
    tmp = J+LAM_1/mu1;
    thr = lambda2/mu1;
    
    if l0norm  % hard thresholding (l0 norm for sparsity)
        thr = sqrt(2*thr);
        if double_thr
            thr = thr*2;
        end
        C2=tmp.*((sign(abs(tmp)-thr)+1)/2);
        
    elseif gmc
        C2 = arrayfun(@(y) firm_thresh(y, thr, thr/gamma), tmp);
    end
    C2 = C2 - diag(diag(C2));
    
    [cm indmax] = maxk(abs(C2),d,1);
    tmp = zeros(N,N);
    
    for n=1:N
        tmp(indmax(:,n),n)=C2(indmax(:,n),n);
    end
    C2 = tmp; clear tmp;
    
    C = lambda1*C1 + lambda2*C2;
    %C = C - diag(diag(C));
    
    % compute Lagrangian
    %sparse_term = sum(C2(:)~=0);
    %low_rank_term = sum(sign(sig));
    
    %L(iter) = 1/2*norm(X-X*J,'fro')^2 + lambda1*low_rank_term...
    %   + lambda2*sparse_term;
    
    %if iter>1
    %   Lgrad = abs(L(iter)-L(iter-1));
    %else
    %   Lgrad = L(iter);
    %end
    
    % Update of Lagrange multipliers
    LAM_1 = LAM_1 + mu1*(J - C);
    
    mu1 = min(rho*mu1, max_mu1);
    
    if rho~=1 || iter==1
        Jf = inv(XT + mu1*eye(N));
    end
    
    %L_grad_vect = [L_grad_vect ; Lgrad];
    
    err1 = max(max(abs(J-C)));
    err2 = max(max(abs(J-J_prev)));
    
    
    err1_all = [err1_all ; err1];
    err2_all = [err2_all ; err2];
    
    if err1<err_thr && err2<err_thr
        not_converged=0;
    end
    
    % check convergence
    if iter>=iter_max
        not_converged = 0;
    end
    
    %if Lgrad < err_thr
    %  not_converged = 0;
    %end
    
    iter = iter+1;
    
end
%iter

%plot(1:length(err1_all), err1_all, '.-');
%figure;
%plot(1:length(err2_all), err2_all, '.-');
%plot(1:length(L_grad_vect), L_grad_vect, '.-');
%L_grad_vect

C = normc(C);
%C=C1;

error = norm(X-X*J)/norm(X);

