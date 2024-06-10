function [labels_est_X, labels_est_Xwp] = run_SC_algorithm(X,Xwp,labels,paras)

% Runs selected subspace clustering algorithm on data in amninet space (X) 
% and in wavelet packets subband (Xwp)
%
%  Inputs:
%

algorithm = paras.algorithm;
nc=max(labels);   % number of clusters

if strcmp (algorithm,'SSC')
    fprintf('Running SSC ...\n');

    outlier = paras.outlierAmbient; affine = paras.affineAmbient;  r= 0; rho = 2.0;
    ipd_ambient_domain=paras.ipd_ambient_domain;
    alpha_X = paras.alphaAmbient;
    [Z_x,labels_est_X] = SSC(normc(X),r,affine,alpha_X,outlier,rho,labels);
    if ipd_ambient_domain
        d_x = paras.dimSubspaceAmbient;
        C_sym = BuildAdjacency_cut(Z_x,d_x);
        labels_est_X = SpectralClusteringL(C_sym,nc);
    end

    outlier = paras.outlierWP; affine = paras.affineWP;  r= 0; rho = 2.0;
    ipd_wp_domain=paras.ipd_wp_domain;
    alpha_wp = paras.alphaWP;
    [Z_wp,labels_est_Xwp] = SSC(normc(Xwp),r,affine,alpha_wp,outlier,rho,labels);
    if ipd_wp_domain
        d_wp =paras.dimSubspaceWP;
        C_sym = BuildAdjacency_cut(Z_wp,d_wp);
        labels_est_Xwp = SpectralClusteringL(C_sym,nc);
    end
elseif strcmp (algorithm,'S0L0_LRSSC')
    fprintf('Running LRSSC l0+S0 with proximal average..\n');

    ipd_ambient_domain=paras.ipd_ambient_domain;
    ipd_wp_domain=paras.ipd_wp_domain; 

    alpha_x = paras.alphaAmbient; lambda_x = paras.lambdaAmbient;
    options = struct('lambda',lambda_x,'alpha',alpha_x, 'err_thr',1e-4,'iter_max',100, 'double_thr', false, ...
        'l0norm', true);

    [Z_x, error] = ADMM_LRSSC_prox_avg(normc(X),options);  

    if ipd_ambient_domain
        d_x = paras.dimSubspaceAmbient;
        C_sym = BuildAdjacency_cut(abs(Z_x),d_x);
        labels_est_X = SpectralClusteringL(C_sym,nc);
    else
        labels_est_X = SpectralClusteringL(abs(Z_x)+abs(Z_x'),nc);
    end

    alpha_wp = paras.alphaWP; lambda_wp = paras.lambdaWP;
    options = struct('lambda',lambda_wp,'alpha',alpha_wp, 'err_thr',1e-4,'iter_max',100, 'double_thr', false, ...
        'l0norm', true);

    [Z_wp, error] = ADMM_LRSSC_prox_avg(normc(Xwp),options);
    if ipd_wp_domain
        d_wp = paras.dimSubspaceWP;
        C_sym = BuildAdjacency_cut(abs(Z_wp),d_wp);
        labels_est_Xwp = SpectralClusteringL(C_sym,nc);
    else
        labels_est_Xwp = SpectralClusteringL(abs(Z_wp)+abs(Z_wp'),nc);
    end
elseif strcmp(algorithm,'GMC_LRSSC')
    fprintf('Running GMC_LRSSC ...\n');

    alpha_x = paras.alphaAmbient; lambda_x = paras.lambdaAmbient;
    gamma_x = paras.gammaAmbient;

    ipd_ambient_domain=paras.ipd_ambient_domain;
    ipd_wp_domain=paras.ipd_wp_domain;

    options = struct('lambda',lambda_x,'alpha',alpha_x,'rank_est',0.6,'gamma',gamma_x,...
        'err_thr',1e-4,'iter_max',100, 'affine',false,...
        'l1_nucl_norm',false,'l0norm',false,'elra',false, 'gmc',true);
    %
    [Z_x, error] = ADMM_LRSSC(normc(X),options);
    if ipd_ambient_domain
        d_x = paras.dimSubspaceAmbient;
        C_sym = BuildAdjacency_cut(abs(Z_x),d_x);
        labels_est_X = SpectralClusteringL(C_sym,nc);
    else
        labels_est_X = SpectralClusteringL(abs(Z_x)+abs(Z_x'),nc);
    end

    alpha_wp = paras.alphaWP; lambda_wp = paras.lambdaWP;
    gamma_wp = paras.gammaWP;
    options = struct('lambda',lambda_wp,'alpha',alpha_wp,'rank_est',0.6,'gamma',gamma_wp,...
        'err_thr',1e-4,'iter_max',100, 'affine',false,...
        'l1_nucl_norm',false,'l0norm',false,'elra',false, 'gmc',true);

    [Z_wp, error] = ADMM_LRSSC(normc(Xwp),options);
    if ipd_wp_domain
        d_wp = paras.dimSubspaceWP;
        C_sym = BuildAdjacency_cut(abs(Z_wp),d_wp);
        labels_est_Xwp = SpectralClusteringL(C_sym,nc);
    else
        labels_est_Xwp = SpectralClusteringL(abs(Z_wp)+abs(Z_wp'),nc);
    end
elseif strcmp(algorithm,'NSN')
    fprintf('Running NSN+Spectral..\n');

    ipd_ambient_domain=paras.ipd_ambient_domain;
    ipd_wp_domain=paras.ipd_wp_domain;

    nbh_x = paras.kAmbient;
    dmax_x = paras.dmaxAmbient;
    Z_x = NSN(X,nbh_x,dmax_x,1e-4);
    if ipd_ambient_domain
        d_x = paras.dimSubspaceAmbient;
        C_sym = BuildAdjacency_cut(abs(Z_x),d_x);
        labels_est_X = SpectralClusteringL(C_sym,nc);
    else
        labels_est_X = SpectralClusteringL(abs(Z_x)+abs(Z_x'),nc);
    end

    nbh_wp = paras.kWP;
    dmax_wp = paras.dmaxWP;    
    Z_wp = NSN(Xwp,nbh_wp,dmax_wp,1e-4);
    if ipd_wp_domain
        d_wp = paras.dimSubspaceWP;
        C_sym = BuildAdjacency_cut(abs(Z_wp),d_wp);
        labels_est_Xwp = SpectralClusteringL(C_sym,nc);
    else
        labels_est_Xwp = SpectralClusteringL(abs(Z_wp)+abs(Z_wp'),nc);
    end
elseif strcmp(algorithm,'RTSC')
    fprintf('Running TSC+Spectral..\n');

    ipd_ambient_domain=paras.ipd_ambient_domain;
    ipd_wp_domain=paras.ipd_wp_domain;

    Nc = size(X,2)/nc;% number of data per cluster

    q=paras.qAmbient;
    q_x = max(q,ceil(Nc/20));
    [labels_est_X,Z_x] = TSC(normc(X),q_x,nc);
    if ipd_ambient_domain
        d_x = paras.dimSubspaceAmbient;
        C_sym = BuildAdjacency_cut(abs(Z_x),d_x);
        labels_est_X = SpectralClusteringL(C_sym,nc);
    end

    q=paras.qWP;
    q_wp = max(q,ceil(Nc/20));
    [labels_est_Xwp,Z_wp] = TSC(normc(Xwp),q_wp,nc);
    if ipd_wp_domain
        d_wp = paras.dimSubspaceWP;
        C_sym = BuildAdjacency_cut(abs(Z_wp),d_wp);
        labels_est_Xwp = SpectralClusteringL(C_sym,nc);
    end
elseif strcmp(algorithm,'LRR')
    fprintf('Running LRR..\n');

    ipd_ambient_domain=paras.ipd_ambient_domain;
    ipd_wp_domain=paras.ipd_wp_domain;

    lambda_x = paras.lambdaAmbient;
    Z_x = solve_lrr(normc(X),lambda_x);
    if ipd_ambient_domain
        d_x = paras.dimSubspaceAmbient;
        Z_x = BuildAdjacency_cut(Z_x,d_x); % keep largest d_x coefficients
    end

    % post processing
    [U,S,V] = svd(Z_x,'econ');
    S = diag(S);
    r = sum(S>1e-4*S(1));
    U = U(:,1:r);S = S(1:r);
    U = U*diag(sqrt(S));
    U = normr(U);
    U = U./repmat(sqrt(sum(U.^2,2)),1,size(U,2));
    LL = (U*U').^4;

    % spectral clustering
    D = diag(1./sqrt(sum(LL,2)));
    LL = D*LL*D;
    [U,S,V] = svd(LL);
    V = U(:,1:nc);
    V = D*V;
    labels_est_X = kmeans(V,nc,'emptyaction','singleton','replicates',10,'display','off');

    % The following scripts are copied from the LRR source code.
    lambda_wp = paras.lambdaWP;
    % run lrr
    Z_wp = solve_lrr(normc(Xwp),lambda_wp);
    if ipd_wp_domain
        d_wp = paras.dimSubspaceWP;
        Z_wp = BuildAdjacency_cut(Z_wp,d_wp); % keep largest d_wp coefficients
    end

    % post processing
    [U,S,V] = svd(Z_wp,'econ');
    S = diag(S);
    r = sum(S>1e-4*S(1));
    U = U(:,1:r);S = S(1:r);
    U = U*diag(sqrt(S));
    U = normr(U);
    U = U./repmat(sqrt(sum(U.^2,2)),1,size(U,2));
    LL = (U*U').^4;

    % spectral clustering
    D = diag(1./sqrt(sum(LL,2)));
    LL = D*LL*D;
    [U,S,V] = svd(LL);
    V = U(:,1:nc);
    V = D*V;
    labels_est_Xwp = kmeans(V,nc,'emptyaction','singleton','replicates',10,'display','off');
elseif strcmp(algorithm,'S_1o2_LRR')
    fprintf('Running LRR_S_1/2..\n');

    ipd_ambient_domain=paras.ipd_ambient_domain;
    ipd_wp_domain=paras.ipd_wp_domain;

    lambda_x = paras.lambdaAmbient;
    Z_x = Run_SpNM('spdual12lrr',normc(X),lambda_x);
    if ipd_ambient_domain
        d_x = paras.dimSubspaceAmbient;
        C_sym = BuildAdjacency_cut(Z_x,d_x);
        labels_est_X = SpectralClusteringL(C_sym,nc);
    else
        labels_est_X = SpectralClusteringL(abs(Z_x)+abs(Z_x'),nc);
    end

    lambda_wp = paras.lambdaWP;
    Z_wp = Run_SpNM('spdual12lrr',normc(Xwp),lambda_wp);
    if ipd_wp_domain
        d_wp = paras.dimSubspaceWP;
        C_sym = BuildAdjacency_cut(Z_wp,d_wp);
        labels_est_Xwp = SpectralClusteringL(C_sym,nc);
    else
        labels_est_Xwp = SpectralClusteringL(abs(Z_wp)+abs(Z_wp'),nc);
    end
elseif strcmp(algorithm,'S_2o3_LRR')
    fprintf('Running LRR_S_2/3..\n');
%     
    ipd_ambient_domain=paras.ipd_ambient_domain;
    ipd_wp_domain=paras.ipd_wp_domain;

   lambda_x = paras.lambdaAmbient;
   Z_x = Run_SpNM('spdual23lrr',normc(X),lambda_x); 
    if ipd_ambient_domain
        d_x = paras.dimSubspaceAmbient;
        C_sym = BuildAdjacency_cut(Z_x,d_x);
        labels_est_X = SpectralClusteringL(C_sym,nc);
    else
        labels_est_X = SpectralClusteringL(abs(Z_x)+abs(Z_x'),nc);
    end
 
    lambda_wp = paras.lambdaWP;
    Z_wp = Run_SpNM('spdual23lrr',normc(Xwp),lambda_wp);
    if ipd_wp_domain
        d_wp = paras.dimSubspaceWP;
        C_sym = BuildAdjacency_cut(Z_wp,d_wp);
        labels_est_Xwp = SpectralClusteringL(C_sym,nc);
    else
        labels_est_Xwp = SpectralClusteringL(abs(Z_wp)+abs(Z_wp'),nc);
    end
end

end
