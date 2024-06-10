%% Wavelet packetb MERA tensor network subspace clustering
% I. Kopriva, D. Seršić 2024-05

%%
clear
close all

% Parameters of 2D wavelet packet decomposition
numberOfLevels = 1;
wavelet = 'haar';  % 'fk4' or 'haar'

% Set path to all subfolders
addpath(genpath('.'));

% Number of evaluations
numit =10;

%% Load the data from the chosen dataset

% Please uncomment the dataset that you want to use and comment the other ones
% dataName = 'YaleBCrop025';
% dataName = 'MNIST';
% dataName = 'USPS';
dataName = 'ORL';
% dataName = 'COIL20'; 
% dataName = 'COIL100'; 

%% prepare data 
if strcmp(dataName,'YaleBCrop025')
   % Dataset dependent parameters
   i1 = 48; i2 = 42; % image size
   dimSubspace=9;   % subspaces dimension
   numIn = 43;  % number of in-sample data
   numOut = 21; % number of out-of-sample data
   ni = 64; % number of images per group
   nc = 38;  % number of groups

  % parameters of MERA network
  rX= [43 38 43 38 5];
  R_opt=10;  % Level 1
  lambda_opt = 1e-3;  

   % postprocessing of coefficient representation matrix
   ipd_postprocessing = false;

   load YaleBCrop025.mat;
   [i1, i2, ni, nc] = size(I); % rectangular image size: i1 x i2, number of images per person: ni, number of persons: nc
   clear I Ind s ns
   
   N = nc*ni; % number of samples
   X = zeros(i1*i2, N); labels = zeros(N,1); % allocation of space
    
   ns = 0; % sample number counter
   for i=1:nc % person
       for j=1:ni % face image
           ns = ns + 1; % sample index
           X(:,ns) = Y(:,j,i); % sample (columns of X represent vectorized data of rectangular images)
           labels(ns,1) = i;    % to be used for oracle based validation 
       end
   end
   Y = X;
   clear X i j
   
elseif strcmp(dataName,'MNIST')
   % Dataset dependent parameters
   i1 = 28; i2 = 28; % image size
   dimSubspace=12;   % subspaces dimension
   numIn = 50;  % number of in-sample data
   numOut = 50; % number of out-of-sample data
   nc = 10;   % number of groups

  % parameters of MERA network
  rX= [25 20 20 25 5];
  R_opt=3;  % Level 1
  lambda_opt = 1e-4;

   % postprocessing of coefficient representation matrix
   ipd_postprocessing = false;
    
   images = loadMNISTImages('t10k-images.idx3-ubyte'); % columns of X represent vectorized data of squared images
   % i1 = 28; i2 = 28; N = 10000; % 10000 images of ten digits (each image 28x28 pixels)
   
   labels = loadMNISTLabels('t10k-labels.idx1-ubyte'); % to be used for oracle based validation
   % nc = 10; % ten hand-written digits
   [labelssorted,IX] = sort(labels);
   Y = images(:,IX);  
   labels = labelssorted + 1; 
   clear labelssorted IX images  

elseif strcmp(dataName,'USPS')
   % Dataset dependant parameters
   i1 = 16; i2 = 16; % image size
   dimSubspace=12;   % subspaces dimension
   numIn = 50;  % number of in-sample data
   numOut = 50; % number of out-of-sample data
   nc = 10;    % number of groups

  % parameters of MERA network
   rX= [25 20 20 25 5];
  R_opt=6;  % Level 1
  lambda_opt = 1e-2;
     
   % postprocessing of coefficient representation matrix
   ipd_postprocessing = false;

   data = load('usps');
   X = data(:,2:end)'; % columns of X represent vectorized data of squared images
   % i1 = 16; i2 = 16; N = 7291; nc = 10; % 1000 images of each of ten digits (each image 16x16 pixels)
    
   labels = data(:,1)-1;  % to be used for oracle based validation
   % nc = 10; % ten hand-written digits
   [labelssorted,IX] = sort(labels);
   Y = X(:,IX);  
   labels = labelssorted + 1; 
   clear data X labelssorted IX

elseif strcmp(dataName,'ORL')  
   % Dataset dependant parameters
   i1 = 32; i2 = 32; % each image 32x32 pixels
   dimSubspace=5;   % subspaces dimension
   numIn = 7;  % number of in-sample data
   numOut = 3; % number of out-of-sample data
   nc = 40; % number of groups

   % parameters of MERA network
   rX= [20 14 20 14 5]; % level 1
   R_opt=10;  % Level 1
   lambda_opt = 1e-9;  

   % postprocessing of coefficient representation matrix
   ipd_postprocessing = true;
  
   data = load('ORL_32x32.mat');  
   Y = data.fea'; % columns of X represent vectorized data of squared images
   % i1 = 32; i2 = 32; N = 400; nc = 40; % 400 face images of 40 persons (each image 32x32 pixels)
    
   labels=data.gnd;   % to be used for oracle based validation
   clear data  

elseif strcmp(dataName,'COIL20')
   % Dataset dependant parameters
   i1=32; i2=32; % image dimensions
   dimSubspace=10;   % subspaces dimension
   numIn = 22;  % number of in-sample data
   numOut = 22; % number of out-of-sample data
   nc = 20; % number of groups

   % parameters of MERA network
  rX= [22 20 22 20 5];  % ??????????????????????????????????????????
  R_opt=13;  % Level 1
  lambda_opt = 0.1;

   % postprocessing of coefficient representation matrix
   ipd_postprocessing = true;  
  
   load COIL20.mat
   Y=transpose(fea); % columns of X represent vectorized data of squared images
   % i1=32; i2=32; N=1440; nc=20; % 1440 images of 20 objects (72 images per object) (each image is 32x32 pixels)
   clear fea;
    
   labels=gnd;   % to be used for oracle based validation
   % nc = 20; % twenty objects images
   clear gnd  

elseif strcmp(dataName,'COIL100')
   % Dataset dependant parameters
   i1=32; i2=32; % image dimensions
   dimSubspace=10;   % subspaces dimension
   numIn = 22;  % number of in-sample data
   numOut = 22; % number of out-of-sample data   
   nc = 100; % number of groups

   % parameters of MERA network
   rX= [50 44 50 44 5];;  % ??????????????????????????????????????????
   R_opt=17;  % Level 1
   lambda_opt = 0.1;
   
   % postprocessing of coefficient representation matrix
   ipd_postprocessing = true;

   load COIL100.mat
   Y=double(fea.'); % columns of X represent vectorized data of squared images
   % i1=32; i2=32; N=7200; nc=100; % 7200 images of 100 objects (72 images per object) (each image is 32x32 pixels)
   clear fea;

   labels=gnd;    % to be used for oracle based validation
   % nc = 100; % one hundred objects images
   clear gnd
end

for it=1:numit
    fprintf('Iter %d\n',it);

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

 %% STEP 2: MERA
    XMV{1} = X_in;
    for sub_ind=1:K
        Xwp_in = squeeze(XWP(:,sub_ind,:));  % to be used by selected SC algorithm
        XMV{sub_ind+1} = Xwp_in;
    end

    nV=length(XMV);
    for v=1:nV
        [XMV{v}]=NormalizeData(XMV{v});
    end

    % parameters of MERA network
    paras_mera.rX = rX;
    paras_mera.R{1}=[R_opt,R_opt];
    paras_mera.lambda=lambda_opt;

    [S,mera]=MERA_MSC(XMV,paras_mera);
    if ipd_postprocessing
       S = BuildAdjacency_cut(S,dimSubspace);
    end

    labels_est(1,:) = SpectralClusteringL(S,nc); 
    ACC_in(it)  = 1 - computeCE(labels_est,labels_in)
    NMI_in(it) = compute_nmi(labels_in,labels_est)
    Rand_in(it) = RandIndex(labels_in,labels_est)
    Fscore_in(it) = compute_f(labels_in,labels_est)
    Purity_in(it) = purFuc(labels_in,labels_est) 

    % Out-of-sample clustering
    for k=1:5
        [B_x(:,:,k), begB_x(1,:,k), endB_x(1,:,k), mu_X(:,:,k)] = bases_estimation(XMV{k},labels_est,dimSubspace);
    end

    % Clustering out-of-sample data
    XWP_out = zeros(rows*cols, K+1, size(X_out,2)); % Space for the result
    XWP_out(:,1,:) = X_out;

    for ns=1:size(X_out,2)
        xn = reshape(X_out(:,ns), i1, i2); % Back to 2D for WP decomposition

        %%%% WP TRANSFORM
        [SWPC, ~] = swpdec2a(xn, numberOfLevels, wavelet);
        for k = 1:K % subbands
            XWP_out(:, k+1, ns) = reshape(SWPC(:,:,k), rows*cols, 1); % STORE in format: XWP(vectorized coefficients, subband_index,ns);
        end
    end

    % calculate distances subband-wise
    A0=labels_out;
    N_out = size(X_out,2);
    for k=1:K+1
        Xwp_out = squeeze(XWP_out(:,k,:));
        Xwp_out = normc(Xwp_out);

        for l=1:nc
            X_outm = Xwp_out - mu_X(:,l,k);    % make data zero mean for distance calculation
            BB=B_x(:,begB_x(l):endB_x(l),k);
            Xproj = (BB*BB')*X_outm;
            Dproj = X_outm - Xproj;
            D(l,:) = sqrt(sum(Dproj.^2,1));
        end
        [Dmin(k,:), Ax(:,k)] = min(D);
        clear D
    end

    [dmin amin] = min(Dmin);
    for ii=1:N_out
        AAx(ii) = Ax(ii,amin(ii));
    end

    % Performance on out-of-sample data
    ACC_out(it)  = 1 - computeCE(AAx,A0)
    NMI_out(it) = compute_nmi(A0,AAx)
    Fscore_out(it) = compute_f(A0,AAx)
    Rand_out(it) = RandIndex(A0,AAx)
    Purity_out(it) = purFuc(A0,AAx)

    clear AAx Ax XMV X_in X_out labels_in labels_out 
end

display('Estimated performances:')

display('*********** In-sample data:')
mean_ACC_in=mean(ACC_in)
std_ACC_in=std(ACC_in)
  
mean_NMI_in=mean(NMI_in)
std_NMI_in=std(NMI_in)

mean_Fscore_in=mean(Fscore_in)
std_Fscore_in=std(Fscore_in)

mean_Rand_in=mean(Rand_in)
std_Rand_in=std(Rand_in)

mean_purity_in=mean(Purity_in)
std_purity_in=std(Purity_in)


display('*********** Out_of sample data:')
mean_ACC_out= mean(ACC_out)
std_ACC_out = std(ACC_out)

mean_NMI_out = mean(NMI_out)
std_NMI_out = std(NMI_out)

mean_Fscore_out=mean(Fscore_out)
std_Fscore_out=std(Fscore_out)

mean_Rand_out=mean(Rand_out)
std_Rand_out=std(Rand_out)

mean_purity_out=mean(Purity_out)
std_purity_out=std(Purity_out)
