function [paras_data,paras_SC] = params_data_and_algorithms(dataName,algorithm)
%
% (c) Ivica Kopriva, May 2024
%
% This function returns parameters for chosen dataset in structure paras_data,
% and parameters for chosen algorithm in structure paras_SC.
%

if strcmp(dataName,'YaleBCrop025')
    % Dataset dependent parameters
    paras_data.i1 = 48; paras_data.i2 = 42; % image size
    paras_data.dimSubspace=9; % % subspaces dimension
    paras_data.numIn = 43; % number of in-sample data
    paras_data.numOut = 21; % number of out-of-sample data
    paras_data.nc = 38; % number of groups

    % Set algorithm dependent (hyper)parameters
    paras_SC.algorithm = algorithm;
    if strcmp(algorithm,'SSC')
        paras_SC.sub_ind = 4; % D subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 8;
        paras_SC.outlierAmbient = true;
        paras_SC.affineAmbient = false;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 3;
        paras_SC.outlierWP = true;
        paras_SC.affineWP = false;
        paras_SC.alphaAmbient = 15;
        paras_SC.alphaWP = 10;
    elseif strcmp(algorithm,'S0L0_LRSSC')
        paras_SC.sub_ind = 17; % DA subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 6;
        paras_SC.lambdaAmbient = 0;
        paras_SC.lambdaWP = 0.5;
        paras_SC.alphaAmbient = 3;
        paras_SC.alphaWP = 9;
    elseif strcmp(algorithm,'GMC_LRSSC')
        paras_SC.sub_ind = 18; % DH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 8;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 9;
        paras_SC.lambdaAmbient = 29;
        paras_SC.lambdaWP = 0.1;
        paras_SC.alphaAmbient = 2.5;
        paras_SC.alphaWP = 0.6;
        paras_SC.gammaAmbient = 1.0;
        paras_SC.gammaWP = 0.6;
    elseif strcmp(algorithm,'NSN')
        paras_SC.sub_ind = 18; % DH subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.ipd_wp_domain = false;
        paras_SC.kAmbient = 32;
        paras_SC.kWP = 42;
        paras_SC.dmaxAmbient = 32;
        paras_SC.dmaxWP = 42;
    elseif strcmp(algorithm,'RTSC')
        paras_SC.sub_ind = 18; % DH subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.ipd_wp_domain = false;
        paras_SC.qAmbient = 2;
        paras_SC.qWP = 7;
    elseif strcmp (algorithm,'LRR')
        paras_SC.sub_ind = 4; % D subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 9;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 9;
        paras_SC.lambdaAmbient = 2;
        paras_SC.lambdaWP = 0.7;
    elseif strcmp(algorithm,'S_1o2_LRR')
        paras_SC.sub_ind = 17; % DA subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 9;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 8;
        paras_SC.lambdaAmbient = 6;
        paras_SC.lambdaWP = 0.6;
    elseif strcmp(algorithm,'S_2o3_LRR')
        paras_SC.sub_ind = 4; % D subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 9;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 7;
        paras_SC.lambdaAmbient = 4;
        paras_SC.lambdaWP = 0.8;
    end

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
    paras_data.X=X;
    paras_data.labels=labels;

elseif strcmp(dataName,'MNIST')
    % Dataset dependant parameters
    paras_data.i1 = 28; paras_data.i2 = 28; % image size
    paras_data.dimSubspace=12; % % subspaces dimension
    paras_data.numIn = 50; % number of in-sample data
    paras_data.numOut = 50; % number of out-of-sample data
    paras_data.nc = 10; % number of groups

    % Set algorithm dependent (hyper)parameters
    paras_SC.algorithm = algorithm;
    if strcmp(algorithm,'SSC')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 10;
        paras_SC.outlierAmbient = true;
        paras_SC.affineAmbient = true;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 6;
        paras_SC.outlierWP = true;
        paras_SC.affineWP = false;
        paras_SC.alphaAmbient = 6;
        paras_SC.alphaWP = 5;
    elseif strcmp(algorithm,'S0L0_LRSSC')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 9;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 10;
        paras_SC.lambdaAmbient = 0.3;
        paras_SC.lambdaWP = 0.35;
        paras_SC.alphaAmbient = 24;
        paras_SC.alphaWP = 9;
    elseif strcmp(algorithm,'GMC_LRSSC')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 12;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 12;
        paras_SC.lambdaAmbient = 0.1;
        paras_SC.lambdaWP = 0.1;
        paras_SC.alphaAmbient = 46;
        paras_SC.alphaWP = 50;
        paras_SC.gammaAmbient = 0.7;
        paras_SC.gammaWP = 1.0;
    elseif strcmp(algorithm,'NSN')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.ipd_wp_domain = false;
        paras_SC.kAmbient = 36;
        paras_SC.kWP = 31;
        paras_SC.dmaxAmbient = 12;
        paras_SC.dmaxWP = 12;
    elseif strcmp(algorithm,'RTSC')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 8;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 6;
        paras_SC.qAmbient = 18;
        paras_SC.qWP = 6;
    elseif strcmp (algorithm,'LRR')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 11;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 11;
        paras_SC.lambdaAmbient = 0.4;
        paras_SC.lambdaWP = 0.3;
    elseif strcmp(algorithm,'S_1o2_LRR')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 9;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 12;
        paras_SC.lambdaAmbient = 0.25;
        paras_SC.lambdaWP = 0.15;
    elseif strcmp(algorithm,'S_2o3_LRR')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 9;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 10;
        paras_SC.lambdaAmbient = 0.15;
        paras_SC.lambdaWP = 0.2;
    end

    X = loadMNISTImages('t10k-images.idx3-ubyte'); % columns of X represent vectorized data of squared images
    labels = loadMNISTLabels('t10k-labels.idx1-ubyte'); % to be used for oracle based validation
    [labelssorted,IX] = sort(labels);

    labels = labelssorted + 1;
    paras_data.X=X(:,IX);
    paras_data.labels=labels;

elseif strcmp(dataName,'USPS')
    % Dataset dependant parameters
    paras_data.i1 = 16; paras_data.i2 = 16; % image size
    paras_data.dimSubspace=12; % % subspaces dimension
    paras_data.numIn = 50; % number of in-sample data
    paras_data.numOut = 50; % number of out-of-sample data
    paras_data.nc = 10; % number of groups

    % Set algorithm dependent (hyper)parameters
    paras_SC.algorithm = algorithm;
    if strcmp(algorithm,'SSC')
        paras_SC.sub_ind = 1; % A subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 11;
        paras_SC.outlierAmbient = false;
        paras_SC.affineAmbient = false;
        paras_SC.ipd_wp_domain = false;
        paras_SC.outlierWP = false;
        paras_SC.affineWP = false;
        paras_SC.alphaAmbient = 3;
        paras_SC.alphaWP = 3;
    elseif strcmp(algorithm,'S0L0_LRSSC')
        paras_SC.sub_ind = 1; % A subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.ipd_wp_domain = false;
        paras_SC.lambdaAmbient = 0.35;
        paras_SC.lambdaWP = 0.35;
        paras_SC.alphaAmbient = 9;
        paras_SC.alphaWP = 8;
    elseif strcmp(algorithm,'GMC_LRSSC')
        paras_SC.sub_ind = 1; % A subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 12;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 12;
        paras_SC.lambdaAmbient = 0.1;
        paras_SC.lambdaWP = 0.1;
        paras_SC.alphaAmbient = 13;
        paras_SC.alphaWP = 11;
        paras_SC.gammaAmbient = 0.8;
        paras_SC.gammaWP = 0.9;
    elseif strcmp(algorithm,'NSN')
        paras_SC.sub_ind = 1; % A subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.ipd_wp_domain = false;
        paras_SC.kAmbient = 60;
        paras_SC.kWP = 54;
        paras_SC.dmaxAmbient = 11;
        paras_SC.dmaxWP = 11;
    elseif strcmp(algorithm,'RTSC')
        paras_SC.sub_ind = 1; % A subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 12;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 12;
        paras_SC.qAmbient = 23;
        paras_SC.qWP = 23;
    elseif strcmp (algorithm,'LRR')
        paras_SC.sub_ind = 1; % A subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 12;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 12;
        paras_SC.lambdaAmbient = 0.2;
        paras_SC.lambdaWP = 0.2;
    elseif strcmp(algorithm,'S_1o2_LRR')
        paras_SC.sub_ind = 1; % A subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 12;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 12;
        paras_SC.lambdaAmbient = 0.16;
        paras_SC.lambdaWP = 0.16;
    elseif strcmp(algorithm,'S_2o3_LRR')
        paras_SC.sub_ind = 1; % A subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 14;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 14;
        paras_SC.lambdaAmbient = 0.15;
        paras_SC.lambdaWP = 0.15;
    end

    data = load('usps');
    Y = data(:,2:end)'; % columns of Y represent vectorized data of squared images
    labels = data(:,1);
    clear data

    paras_data.X=Y;
    paras_data.labels=labels;

elseif strcmp(dataName,'ORL')
    % Dataset dependent parameters
    paras_data.i1 = 32; paras_data.i2 = 32; % image size
    paras_data.dimSubspace=5; % % subspaces dimension
    paras_data.numIn = 7; % number of in-sample data
    paras_data.numOut = 3; % number of out-of-sample data
    paras_data.nc = 40; % number of groups

    % Set algorithm dependent (hyper)parameters
    paras_SC.algorithm = algorithm;
    if strcmp(algorithm,'SSC')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 7;
        paras_SC.outlierAmbient = false;
        paras_SC.affineAmbient = true;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 4;
        paras_SC.outlierWP = true;
        paras_SC.affineWP = true;
        paras_SC.alphaAmbient = 19;
        paras_SC.alphaWP = 14;
    elseif strcmp(algorithm,'S0L0_LRSSC')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 5;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 5;
        paras_SC.lambdaAmbient = 0.3;
        paras_SC.lambdaWP = 0.4;
        paras_SC.alphaAmbient = 3;
        paras_SC.alphaWP = 12;
    elseif strcmp(algorithm,'GMC_LRSSC')
        paras_SC.sub_ind = 1; % A subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 5;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 4;
        paras_SC.lambdaAmbient = 1.4;
        paras_SC.lambdaWP = 2.3;
        paras_SC.alphaAmbient = 0.8;
        paras_SC.alphaWP = 1.0;
        paras_SC.gammaAmbient = 0.1;
        paras_SC.gammaWP = 0.5;
    elseif strcmp(algorithm,'NSN')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.ipd_wp_domain = false;
        paras_SC.kAmbient = 6;
        paras_SC.kWP = 3;
        paras_SC.dmaxAmbient = 9;
        paras_SC.dmaxWP = 8;
    elseif strcmp(algorithm,'RTSC')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.ipd_wp_domain = false;
        paras_SC.qAmbient = 3;
        paras_SC.qWP = 4;
    elseif strcmp (algorithm,'LRR')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 7;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 10;
        paras_SC.lambdaAmbient = 1.5;
        paras_SC.lambdaWP = 1.85;
    elseif strcmp(algorithm,'S_1o2_LRR')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 6;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 6;
        paras_SC.lambdaAmbient = 7.5;
        paras_SC.lambdaWP = 20;
    elseif strcmp(algorithm,'S_2o3_LRR')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 6;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 5;
        paras_SC.lambdaAmbient = 4.5;
        paras_SC.lambdaWP = 0.0;
    end

    data = load('ORL_32x32.mat');
    Y = data.fea'; % columns of X represent vectorized data of squared images
    labels=data.gnd;

    paras_data.X=Y;
    paras_data.labels=labels;

elseif strcmp(dataName,'COIL20')
    % Dataset dependant parameters
    paras_data.i1 = 32; paras_data.i2 = 32; % image size
    paras_data.dimSubspace=10; % % subspaces dimension
    paras_data.numIn = 50; % number of in-sample data
    paras_data.numOut = 22; % number of out-of-sample data
    paras_data.nc = 20; % number of groups

    % Set algorithm dependent (hyper)parameters
    paras_SC.algorithm = algorithm;
    if strcmp(algorithm,'SSC')
        paras_SC.sub_ind = 1; % A subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 8;
        paras_SC.outlierAmbient = false;
        paras_SC.affineAmbient = false;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 8;
        paras_SC.outlierWP = false;
        paras_SC.affineWP = false;
        paras_SC.alphaAmbient = 7;
        paras_SC.alphaWP = 12;
    elseif strcmp(algorithm,'S0L0_LRSSC')
        paras_SC.sub_ind = 9; % HA subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 8;
        paras_SC.ipd_wp_domain =false;
        paras_SC.lambdaAmbient = 0.0;
        paras_SC.lambdaWP = 0.5;
        paras_SC.alphaAmbient = 19;
        paras_SC.alphaWP = 24;
    elseif strcmp(algorithm,'GMC_LRSSC')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 6;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 6;
        paras_SC.lambdaAmbient = 10;
        paras_SC.lambdaWP = 0.1;
        paras_SC.alphaAmbient = 6;
        paras_SC.alphaWP = 61;
        paras_SC.gammaAmbient = 0.7;
        paras_SC.gammaWP = 0.5;
    elseif strcmp(algorithm,'NSN')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = false; 
        paras_SC.ipd_wp_domain = true; 
        paras_SC.dimSubspaceWP = 10;        
        paras_SC.kAmbient = 24;
        paras_SC.kWP = 17;
        paras_SC.dmaxAmbient = 24;
        paras_SC.dmaxWP = 17;
    elseif strcmp(algorithm,'RTSC')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 5;
        paras_SC.ipd_wp_domain = false;  
        paras_SC.dimSubspaceWP = 6;        
        paras_SC.qAmbient = 4;
        paras_SC.qWP = 5;
    elseif strcmp (algorithm,'LRR')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 7;
        paras_SC.ipd_wp_domain = true;   
        paras_SC.dimSubspaceWP = 10;
        paras_SC.lambdaAmbient = 0.4;
        paras_SC.lambdaWP = 0.1;
    elseif strcmp(algorithm,'S_1o2_LRR')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 7;
        paras_SC.ipd_wp_domain = false;       
        paras_SC.lambdaAmbient = 0.6;
        paras_SC.lambdaWP = 0.1;
    elseif strcmp(algorithm,'S_2o3_LRR')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 8;
        paras_SC.ipd_wp_domain = false;
        paras_SC.lambdaAmbient = 0.9;
        paras_SC.lambdaWP = 0.1;
    end

    load COIL20.mat
    X=transpose(fea); % columns of X represent vectorized data of squared images
    % i1=32; i2=32; N=1440; nc=20; % 1440 images of 20 objects (72 images per object) (each image is 32x32 pixels)
    clear fea;

    labels=gnd;   % to be used for oracle based validation
    % nc = 20; % twenty objects images
    clear gnd

    paras_data.X=X;
    paras_data.labels=labels;

elseif strcmp(dataName,'COIL100')
    % Dataset dependant parameters
    paras_data.i1=32; paras_data.i2=32; % image dimensions
    paras_data.dimSubspace=10;   % subspaces dimension
    paras_data.numIn = 50;  % number of in-sample data
    paras_data.numOut = 22; % number of out-of-sample data   
    paras_data.nc = 100; % number of groups

     % Set algorithm dependent (hyper)parameters
    paras_SC.algorithm = algorithm;
    if strcmp(algorithm,'SSC')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 2;
        paras_SC.outlierAmbient = false;
        paras_SC.affineAmbient = false;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 3;
        paras_SC.outlierWP = false;
        paras_SC.affineWP = false;
        paras_SC.alphaAmbient = 20;
        paras_SC.alphaWP = 4;
    elseif strcmp(algorithm,'S0L0_LRSSC')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.ipd_wp_domain =false;
        paras_SC.lambdaAmbient = 0.8;
        paras_SC.lambdaWP = 0.4;
        paras_SC.alphaAmbient = 20;
        paras_SC.alphaWP = 22;
    elseif strcmp(algorithm,'GMC_LRSSC')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.ipd_wp_domain = false;
        paras_SC.lambdaAmbient = 9.5;
        paras_SC.lambdaWP = 0.1;
        paras_SC.alphaAmbient = 6;
        paras_SC.alphaWP = 26;
        paras_SC.gammaAmbient = 0.0;
        paras_SC.gammaWP = 0.9;
    elseif strcmp(algorithm,'NSN')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = false; 
        paras_SC.ipd_wp_domain = false;      
        paras_SC.kAmbient = 42;
        paras_SC.kWP = 31;
        paras_SC.dmaxAmbient = 9;
        paras_SC.dmaxWP = 31;
    elseif strcmp(algorithm,'RTSC')
        paras_SC.sub_ind = 5; % AA subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.ipd_wp_domain = false;         
        paras_SC.qAmbient = 2;
        paras_SC.qWP = 4;
    elseif strcmp (algorithm,'LRR')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 9;
        paras_SC.ipd_wp_domain = true;   
        paras_SC.dimSubspaceWP = 4;
        paras_SC.lambdaAmbient = 0.1;
        paras_SC.lambdaWP = 0.6;
    elseif strcmp(algorithm,'S_1o2_LRR')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = false;
        paras_SC.dimSubspaceAmbient = 4;
        paras_SC.ipd_wp_domain = false;
        paras_SC.dimSubspaceWP = 4;
        paras_SC.lambdaAmbient = 0.3;
        paras_SC.lambdaWP = 0.01;
    elseif strcmp(algorithm,'S_2o3_LRR')
        paras_SC.sub_ind = 6; % AH subband
        paras_SC.ipd_ambient_domain = true;
        paras_SC.dimSubspaceAmbient = 4;
        paras_SC.ipd_wp_domain = true;
        paras_SC.dimSubspaceWP = 4;
        paras_SC.lambdaAmbient = 0.15;
        paras_SC.lambdaWP = 0.01;
    end

    load COIL100.mat
    X=double(fea.'); % columns of X represent vectorized data of squared images
    clear fea;

    labels=gnd;    % to be used for oracle based validation
    clear gnd

    paras_data.X=X;
    paras_data.labels=labels;

end

end