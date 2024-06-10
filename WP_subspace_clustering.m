%% Wavelet packet subspace clustering
% I. Kopriva 2024-05
% D. Sersic  2024-05

%% Initialization
clear
close all

% Set path to all subfolders
addpath(genpath('.'));

%% Dataset independant parameters

% Parameters of 2D wavelet packet decomposition
numberOfLevels = 2;
wavelet = 'haar';

% Number of iterations
nIter = 30;

%% Load the data from the chosen dataset

% Please uncomment the dataset that you want to use and comment the other ones
% dataName = 'YaleBCrop025';
% dataName = 'MNIST';
% dataName = 'USPS';
% dataName = 'ORL';
% dataName = 'COIL20'; 
 dataName = 'COIL100'; % Color images

%% Selection of subspace clustering algorithm

% Please uncomment algorithm you want to use and comment the other ones
% algorithm = 'SSC';  % Sparse subspace clustering
% algorithm = 'S0L0_LRSSC'; % S0-L0 constrained low rank sparse SC
% algorithm = 'GMC_LRSSC';  % GMC constrained low rank sparse SC
% algorithm = 'NSN';  % Nearest subspace neighbor SC 
% algorithm = 'RTSC'; % Robust thresholding SC 
% algorithm = 'LRR';  % Nuclear norm low-rank representation SC 
% algorithm = 'S_1o2_LRR'; % L_1/2 norm LRR SC
 algorithm = 'S_2o3_LRR'; % L_2/3 norm LRR SC

%% STEP 1A: prepare data and select hyperparameters of chosen SC algorithm
[paras_data,paras_SC] = params_data_and_algorithms(dataName,algorithm);

i1 = paras_data.i1; i2 = paras_data.i2; % image size
dimSubspace = paras_data.dimSubspace; % Subspace dimension
numIn = paras_data.numIn; % number of in-sample data
numOut = paras_data.numOut; % number of out-of-sample data
nc = paras_data.nc; % number of groups
Y = paras_data.X;  % data
labels = paras_data.labels; % labels
sub_ind = paras_SC.sub_ind; % data set and SC algorithm dependent subband index

ACC_x_in     = zeros(1, nIter);     ACC_x_out    = zeros(1, nIter);
ACC_wp_in    = zeros(1, nIter);     ACC_wp_out   = zeros(1, nIter);
NMI_x_in     = zeros(1, nIter);     NMI_x_out    = zeros(1, nIter);
NMI_wp_in    = zeros(1, nIter);     NMI_wp_out   = zeros(1, nIter);
Fscore_x_in  = zeros(1, nIter);     Fscore_x_out = zeros(1, nIter);
Fscore_wp_in = zeros(1, nIter);     Fscore_wp_out= zeros(1, nIter);
Rand_x_in    = zeros(1, nIter);     Rand_x_out   = zeros(1, nIter);
Rand_wp_in   = zeros(1, nIter);     Rand_wp_out  = zeros(1, nIter);
Purity_x_in  = zeros(1, nIter);     Purity_x_out = zeros(1, nIter);
Purity_wp_in = zeros(1, nIter);     Purity_wp_out= zeros(1, nIter);

affinity_x   = zeros(1, nIter);     affinity_wp = zeros(1, nIter);

for it = 1:nIter
    fprintf('Iter: %d\n',it);

    %% Generate a problem instance
    % Problem instance is a random split of the chosen dataset into an input set (X_in) and an output set (X_out),
    % as well as the concommitant label sets (label_in, label_out)

    rng('shuffle');

    %% STEP 1A prepare in-sample and out-of-sample random partitions
    % Each category is separately split, to ensure proportional representation
    nIn = 1; nOut = 1;
    for c=1:nc % Through all categories
        ind = (labels == c); % Indices of the chosen category
        Xc = Y(:,ind);       % Samples ...
        numSamples = size(Xc, 2); % Number of samples ...
        ind = randperm(numSamples); % Random permutation of the indices
        X_in(:,    nIn:nIn+numIn-1 ) = Xc(:, ind(1:numIn)); % Data
        X_out(:, nOut:nOut+numOut-1) = Xc(:, ind(numIn+1:numIn+numOut));
        labels_in(  nIn:nIn + numIn-1) = c; % Labels
        labels_out(nOut:nOut+numOut-1) = c;
        nIn  = nIn  + numIn; % Next indices
        nOut = nOut + numOut;
    end
    X_in( :,   nIn:end) = []; % Cut out the surplus of the allocated space
    X_out(:,  nOut:end) = [];
    labels_in(  nIn:end) = [];
    labels_out(nOut:end) = [];

    %% STEP 1B: 2D Wavelet packet decomposition

    % Space allocation
    xn = reshape(X_in(:,1), i1, i2);  % Back to 2D for WP decomposition
    [SWPC, ~] = swpdec2a(xn, numberOfLevels, wavelet); % Dummy decomposition
    [rows, cols, K] = size(SWPC);
    XWP = zeros(rows*cols, K, size(X_in,2)); % Space for the result

    for ns=1:size(X_in,2)
        xn = reshape(X_in(:,ns), i1, i2); % Back to 2D for WP decomposition

        %%%% WP TRANSFORM
        [SWPC, ~] = swpdec2a(xn, numberOfLevels, wavelet);
        for k = 1:K % subbands
            XWP(:, k, ns) = reshape(SWPC(:,:,k), rows*cols, 1); % STORE in format: XWP(vectorized coefficients, subband_index,ns);
        end
    end

    Xwp_in = squeeze(XWP(:,sub_ind,:));  % to be used by the selected SC algorithm

    %% STEP 2: apply selected subspace clustering algorithm(s) to Xwp_sub as well as to the original data X_in
    [labels_est_X(1,:), labels_est_Xwp(1,:)] = run_SC_algorithm(X_in,Xwp_in,labels_in,paras_SC);

    %% Performance on in-sample data
    ACC_x_in(it)  = 1 - computeCE(labels_est_X,labels_in)     
    ACC_wp_in(it)  = 1 - computeCE(labels_est_Xwp,labels_in)
    NMI_x_in(it) = compute_nmi(labels_in,labels_est_X)       
    NMI_wp_in(it) = compute_nmi(labels_in,labels_est_Xwp)
    Fscore_x_in(it) = compute_f(labels_in,labels_est_X);         Fscore_wp_in(it) = compute_f(labels_in,labels_est_Xwp);
    Rand_x_in(it) = RandIndex(labels_in,labels_est_X);           Rand_wp_in(it) = RandIndex(labels_in,labels_est_Xwp);
    Purity_x_in(it) = purFuc(labels_in,labels_est_X);            Purity_wp_in(it) = purFuc(labels_in,labels_est_Xwp);

    %% STEP 3: calculate minimal principal angles between the subspaces and estimate the ratio

    % estimate bases in the input ambient space
    if paras_SC.ipd_ambient_domain
        dimSubspaceA = paras_SC.dimSubspaceAmbient;
    else
        dimSubspaceA = dimSubspace;
    end
    [affinity_x(it), B_x, begB_x, enddB_x, mu_X]  = average_affinity(X_in,labels_est_X,dimSubspaceA);
    %
    
    %       % estimate bases in the WP space
    if paras_SC.ipd_wp_domain
        dimSubspaceWP = paras_SC.dimSubspaceWP;
    else
        dimSubspaceWP = dimSubspace;
    end
    [affinity_wp(it), B_wp, begB_wp, enddB_wp, mu_Xwp]  = average_affinity(Xwp_in,labels_est_Xwp,dimSubspaceWP);

    %% Clustering of out-of-sample data
    XWP_out = zeros(rows*cols, K, size(X_out,2)); % Space for the result

    % Clustering out-of-sample data
    XWP_out = zeros(rows*cols, K, size(X_out,2)); % Space for the result

    for ns=1:size(X_out,2)
        xn = reshape(X_out(:,ns), i1, i2); % Back to 2D for WP decomposition

        %%%% WP TRANSFORM
        [SWPC, ~] = swpdec2a(xn, numberOfLevels, wavelet);
        for k = 1:K % subbands
            XWP_out(:, k, ns) = reshape(SWPC(:,:,k), rows*cols, 1); % STORE in format: XWP(vectorized coefficients, subband_index,ns);
        end
    end

    Xwp_out = squeeze(XWP_out(:,sub_ind,:));  % to be used for out-of-sample clustering
    %
    A0=labels_out;  N_out = size(X_out,2);
    X_out = normc(X_out); Xwp_out = normc(Xwp_out);

    for l=1:nc
        X_outm = X_out - mu_X(:,l);    % make data zero mean for distance calculation
        BB=B_x(:,begB_x(l):enddB_x(l));
        Xproj = (BB*BB')*X_outm;
        Dproj = X_outm - Xproj;
        D(l,:) = sqrt(sum(Dproj.^2,1));
    end
    [~, A_x] = min(D);
    clear D

    for l=1:nc
        Xwp_outm = Xwp_out-mu_Xwp(:,l);  % make data zero mean for distance calculation
        BB=B_wp(:,begB_wp(l):enddB_wp(l));
        Xproj = (BB*BB')*Xwp_outm;
        Dproj = Xwp_outm - Xproj;
        D(l,:) = sqrt(sum(Dproj.^2,1));
    end
    [~, A_wp] = min(D);
    clear D

    % Performance on out-of-sample data with algorithm estiated labels
    ACC_x_out(it)  = 1 - computeCE(A_x,A0); ACC_wp_out(it)  = 1 - computeCE(A_wp,A0);
    NMI_x_out(it) = compute_nmi(A0,A_x); NMI_wp_out(it) = compute_nmi(A0,A_wp);
    Rand_x_out(it) = RandIndex(A0,A_x);  Rand_wp_out(it) = RandIndex(A0,A_wp);
    Fscore_x_out(it) = compute_f(A0,A_x); Fscore_wp_out(it) = compute_f(A0,A_wp);
    Purity_x_out(it) = purFuc(A0,A_x); Purity_wp_out(it) = purFuc(A0,A_wp);

    clear A_x A_wp
    %% Iterations (END LOOP)
end

display('Estimated performances:')

display('*********** In-sample data:')
mean_ACC_x_in=mean(ACC_x_in)
std_ACC_x_in=std(ACC_x_in)
mean_ACC_wp_in=mean(ACC_wp_in)
std_ACC_wp_in=std(ACC_wp_in)
  
mean_NMI_x_in=mean(NMI_x_in)
std_NMI_x_in=std(NMI_x_in)
mean_NMI_wp_in=mean(NMI_wp_in)
std_NMI_wp_in=std(NMI_wp_in)

mean_Fscore_x_in=mean(Fscore_x_in)
std_Fscore_x_in=std(Fscore_x_in)
mean_Fscore_wp_in=mean(Fscore_wp_in)
std_Fscore_wp_in=std(Fscore_wp_in)

mean_Rand_x_in=mean(Rand_x_in)
std_Rand_x_in=std(Rand_x_in)
mean_Rand_wp_in=mean(Rand_wp_in)
std_Rand_wp_in=std(Rand_wp_in)

mean_purity_x_in=mean(Purity_x_in)
std_purity_x_in=std(Purity_x_in)
mean_purity_wp_in=mean(Purity_wp_in)
std_purity_wp_in=std(Purity_wp_in)

mean_affinity_x = mean(affinity_x)
std_affinity_x = std(affinity_x)
mean_affinity_wp = mean(affinity_wp)
std_affinity_wp = std(affinity_wp)


display('*********** Out_of sample data:')
mean_ACC_x_out=mean(ACC_x_out)
std_ACC_x_out=std(ACC_x_out)
mean_ACC_wp_out=mean(ACC_wp_out)
std_ACC_wp_out=std(ACC_wp_out)
  
mean_NMI_x_out=mean(NMI_x_out)
std_NMI_x_out=std(NMI_x_out)
mean_NMI_wp_out=mean(NMI_wp_out)
std_NMI_wp_out=std(NMI_wp_out)

mean_Fscore_x_out=mean(Fscore_x_out)
std_Fscore_x_out=std(Fscore_x_out)
mean_Fscore_wp_out=mean(Fscore_wp_out)
std_Fscore_wp_out=std(Fscore_wp_out)

mean_Rand_x_out=mean(Rand_x_out)
std_Rand_x_out=std(Rand_x_out)
mean_Rand_wp_out=mean(Rand_wp_out)
std_Rand_wp_out=std(Rand_wp_out)

mean_purity_x_out=mean(Purity_x_out)
std_purity_x_out=std(Purity_x_out)
mean_purity_wp_out=mean(Purity_wp_out)
std_purity_wp_out=std(Purity_wp_out)
