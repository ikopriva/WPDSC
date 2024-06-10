clear all; close all; clc;

% Before running this code, SSC and LRR codes from the authors' websites
% should be in the subfolders with the following titles
addpath SSC_ADMM_v1.1
addpath code2
addpath supp_material_RSCT/include
addpath SCC

load YaleBCrop025.mat             % We used the resized raw images provided along with the SSC codes.

Y0 = Y;

p = 2016; n = 64;

n_trial = 100;              % Take 100 random subsets of 38 people for each number of clusters

CE  = zeros(5,100);         % record Clustering Error of 5 algorithms
NSE = zeros(5,100);         % record Neighborhood Selection Error of 5 algorithms
ET  = zeros(5,100);         % record Elapsed Time of 5 algorithms

L = 2
    
cluster_id = zeros(n_trial,L);
for i_trial = 1:n_trial
    cluster_id(i_trial,:) = randsample(38,L);
end
    
for i_trial = 1:n_trial

i_trial
cluster_id(i_trial,:)

%% Choose L random clusters out of 38    
N = n*L;

Y = [];
for i=1:L
    Y = [Y Y0(:,:,cluster_id(i_trial,i))];
end
A0 = reshape(repmat(1:L,n,1),1,N);
 
%% K-means
fprintf('Running K-means..\n'); i_algo = 1;

tic;
[A,~] = kmeans(Y',L,'emptyaction','singleton','replicates',20,'display','off');
ET(i_algo,i_trial)  = toc;
CE(i_algo,i_trial)  = computeCE(A,A0);

%% K-flats
fprintf('Running K-flats..\n'); i_algo = 2;

CE(i_algo,i_trial)  = 1;
for i_replicate = 1:1
    tic;
    A = Kflats(Y,9,L);
    ET(i_algo,i_trial)  = toc;
    CE(i_algo,i_trial)  = min(computeCE(A,A0), CE(i_algo,i_trial));
end

%% SSC
%fprintf('Running SSC..\n'); i_algo = 3;
%r = 0; affine = false; alpha = 20; outlier = true; rho = 1;  % We copied the procedure in "run_SSC_Faces.m" from the authors.
%tic
%[~,Z,A] = SSC(Y,r,affine,alpha,outlier,rho,A0);  % And we modified the original code "SSC.m" so that it returns estimated labels.
%ET(i_algo,i_trial)  = toc;
%CE(i_algo,i_trial)  = computeCE(A,A0);

%% LRR
% fprintf('Running LRR..\n'); i_algo = 4;
% A = lrr_face_seg(normc(Y),A0,L);  % We modified the original code "lrr_face_seg.m"
%                            % from the authors so that we can pass the
%                            % face images and get back the estimated labels.
% CE(i_algo,i_trial) = computeCE(A,A0);

%% SCC
fprintf('Running SCC..\n'); i_algo = 5;
%tic
%[A,~] = scc(Y',9,L);
%ET(i_algo,i_trial)  = toc;
%CE(i_algo,i_trial)  = computeCE(A,A0);

%% SSC-OMP
fprintf('Running SSC-OMP..\n'); i_algo = 6;
tic
Z = OMPSC(normc(Y),5);
A = SpectralClusteringL(abs(Z)+abs(Z)',L);
ET(i_algo,i_trial)  = toc;
CE(i_algo,i_trial)  = computeCE(A,A0);

%% TSC
fprintf('Running TSC..\n'); i_algo = 7;
tic
[A,Z] = TSC(Y,3,L);
ET(i_algo,i_trial)  = toc;
CE(i_algo,i_trial)  = computeCE(A,A0);

%% NSN+Spectral
fprintf('Running NSN+Spectral..\n'); i_algo = 8;
tic
Z = NSN(Y,20,20,1e-4);
A = SpectralClusteringL(Z+Z',L);
ET(i_algo,i_trial)  = toc;
CE(i_algo,i_trial)  = computeCE(A,A0);

%%
if i_trial == 1
    CE(:,1)'
    NSE(:,1)'
    ET(:,1)'
else
    [mean(CE(:,1:i_trial)'); median(CE(:,1:i_trial)'); max(CE(:,1:i_trial)')]
    [mean(NSE(:,1:i_trial)'); median(NSE(:,1:i_trial)'); max(NSE(:,1:i_trial)')]
    [mean(ET(:,1:i_trial)')]
end

pause
end







