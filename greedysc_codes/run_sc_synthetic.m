clear all; close all; clc;

% Before running this code, SSC and LRR codes from the authors' websites
% should be in the subfolders with the following titles
addpath SSC_ADMM_v1.1
addpath code2
addpath supp_material_RSCT/include

n_settings = 25;                  % # experimental settings

p_exp = [ 5  5  5  5  5 10 10 10 10 10 20 20 20 20 20 35 35 35 35 35 50 50 50 50 50];      % ambient dimension
d_exp = [ 3  3  3  3  3  6  6  6  6  6 12 12 12 12 12 21 21 21 21 21 30 30 30 30 30];      % subspace dimension
L_exp = [ 5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5];    % # subspaces
n_exp = [ 2  4  6  8 10  2  4  6  8 10  2  4  6  8 10  2  4  6  8 10  2  4  6  8 10] .* d_exp;

n_algo = 6;                 % # algorithms : SSC, GFS, LRR, TSC, NSN+GSR, NSN+Spectral
n_trial = 1;               

E_exp = zeros(n_algo,n_settings);     % Clustering error
C_exp = zeros(n_algo,n_settings);     % Proportion of successful points

for i_trial = 1:n_trial          % # trials
for i_exp = 1:n_settings
    
    [i_exp i_trial]
    
    %% Parameters
    p = p_exp(i_exp);                       % Ambient dimension
    L = L_exp(i_exp);                       % # subspaces
    d = d_exp(i_exp);                       % subspace dimension
    n = n_exp(i_exp)*ones(1,L);             % # sample points for each subspace
    N = sum(n);                             % Total # of samples
    
    %% True subspace generation   
    D0 = cell(1,L);
    for i=1:L
        D0{i} = orth(randn(p,d));
    end

    %% Data point generation
    A0 = zeros(1,N);           % True labels for sample points
    Y  = zeros(p,N);           % Sample points
    X0 = zeros(d,N);           % Weights for the sample points
    
    IDX = [1 cumsum(n)+1];
    for i=1:L
        X = normc(randn(d,n(i)));   % uniformly random unit vectors
    
        A0(:,IDX(i):IDX(i+1)-1) = i;
         Y(:,IDX(i):IDX(i+1)-1) = D0{i}*X;
        X0(:,IDX(i):IDX(i+1)-1) = X;
    end
    
    %% SSC (l1-min)
    display('SSC..'); i_algo = 1;
    Z = admmLasso_mat_func(Y,false,20); Z_SSC = Z;
    W = abs(Z) + abs(Z)';
    A = SpectralClusteringL(W,L);
    D = EstimateSubspace(Y,A,d,L);
    
    E_exp(i_algo,i_exp) = E_exp(i_algo,i_exp) + computeCE(A,A0) / n_trial;
    C_exp(i_algo,i_exp) = C_exp(i_algo,i_exp) + computeNSE(Z,A0,d) / n_trial;
    
    %% SSC-OMP (OMP)
    display('GFS..'); i_algo = 2;
    Z = OMPSC(Y,0); Z_GFS = Z;
    W = abs(Z) + abs(Z)';
    A = SpectralClusteringL(W,L);
    D = EstimateSubspace(Y,A,d,L);
    
    E_exp(i_algo,i_exp) = E_exp(i_algo,i_exp) + computeCE(A,A0) / n_trial;
    C_exp(i_algo,i_exp) = C_exp(i_algo,i_exp) + computeNSE(Z,A0,d) / n_trial;
      
    %% LRR (nuclear norm minimization)
    display('LRR..'); i_algo = 3;
    Z = solve_lrr(Y,1e5); Z_LRR = Z;
    W = abs(Z) + abs(Z)';
    A = SpectralClusteringL(W,L);
    D = EstimateSubspace(Y,A,d,L);
    
    E_exp(i_algo,i_exp) = E_exp(i_algo,i_exp) + computeCE(A,A0) / n_trial;
    C_exp(i_algo,i_exp) = C_exp(i_algo,i_exp) + computeNSE(Z,A0,d) / n_trial; 
    
    %% TSC (NN)
    display('TSC..'); i_algo = 4; K = d-1;
    [A,Z,~] = TSC(Y,K,L); Z_TSC = Z;
        
    E_exp(i_algo,i_exp) = E_exp(i_algo,i_exp) + computeCE(A,A0) / n_trial;
    C_exp(i_algo,i_exp) = C_exp(i_algo,i_exp) + computeNSE(Z,A0,d) / n_trial;

    %% NSN + Spectral clustering
    display('NSN+Spectral..'); i_algo = 5; K = d-1; kmax = d;
    [Z,~] = NSN(Y,K); Z_NSN = Z;
    W = abs(Z) + abs(Z)';
    A = SpectralClusteringL(W,L);
    D = EstimateSubspace(Y,A,d,L);
        
    E_exp(i_algo,i_exp) = E_exp(i_algo,i_exp) + computeCE(A,A0) / n_trial;
    C_exp(i_algo,i_exp) = C_exp(i_algo,i_exp) + computeNSE(Z,A0,d) / n_trial; 
    
    %% NSN + greedy subspace recovery
    display('NSN+GSR..'); i_algo = 6; K = d-1; kmax = d; epsilon = 1e-4;
    [Z,U] = NSN(Y,K);
    
    D = cell(1,L); I = 1:N; l = 1;
    for l=1:L
        npoints = zeros(1,numel(I)); Y1 = Y(:,I);
        for i=1:numel(I)
            npoints(i) = nnz(sum((Y1-U{I(i)}*(U{I(i)}'*Y1)).^2,1) < epsilon);
        end
        [~,maxindex] = max(npoints);
        
        D{l} = U{I(maxindex)};
        
        J = find(sum((Y1-D{l}*(D{l}'*Y1)).^2,1) < epsilon);
        I = setdiff(I,I(J));
    end
    A = LabelPoints(Y,D);
        
    E_exp(i_algo,i_exp) = E_exp(i_algo,i_exp) + computeCE(A,A0) / n_trial;
    C_exp(i_algo,i_exp) = C_exp(i_algo,i_exp) + computeNSE(Z,A0,d) / n_trial; 
    
    %%
    E_exp
    C_exp
end
end

close all;
n_algo = 6; algo = {'SSC', 'SSC-OMP', 'LRR', 'TSC', 'NSN+Spectral', 'NSN+GSR'};
figure;
set(gcf, 'Position', [100 100 1000 200]);
for i_algo=1:n_algo
    subplot(1,n_algo,i_algo);
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1) pos(2)+.1 pos(3) pos(4)-.2] );
    h = imagesc(flipud(reshape(E_exp(i_algo,:),5,5)'),[0 1]); colormap(flipud(gray));
    set(gca, 'xTick', 1:5);
    set(gca, 'xTickLabel', 2:2:10);
    set(gca, 'yTick', 1:5);
    set(gca, 'yTickLabel', [50 35 20 10 5]);
    title(algo{i_algo});
    if i_algo == 1
        ylabel('Ambient dimension (p)');
    end
    if i_algo == 4
        xlabel('Number of points per dimension for each subspace (n/d)');
    end
end
colorbar;

orient landscape;
set(gcf,'PaperPosition',[.25 3.0 10.5 2.5]);
%print -dpdf test

%

n_algo = 5; algo = {'l1-minimization (SSC)', 'OMP (SSC-OMP)', 'Nuclear norm min. (LRR)', 'Nearest neighbor (TSC)', 'NSN'};
figure;
set(gcf, 'Position', [100 100 1000 200]);
for i_algo=1:n_algo
    subplot(1,n_algo,i_algo);
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1) pos(2)+.1 pos(3) pos(4)-.2] );
    h = imagesc(flipud(reshape(C_exp(i_algo,:),5,5)'),[0 1]); colormap(flipud(gray));
    set(gca, 'xTick', 1:5);
    set(gca, 'xTickLabel', 2:2:10);
    set(gca, 'yTick', 1:5);
    set(gca, 'yTickLabel', [50 35 20 10 5]);
    title(algo{i_algo});
    if i_algo == 1
        ylabel('Ambient dimension (p)');
    end
    if i_algo == 3
        xlabel('Number of points per dimension for each subspace (n/d)');
    end
end
colorbar;

orient landscape;
set(gcf,'PaperPosition',[.25 3.0 10.5 2.5]);
%print -dpdf test2


