clear all; close all; clc;

% Before running this code, SSC and LRR codes from the authors' websites
% should be in the subfolders with the following titles
addpath SSC_ADMM_v1.1
addpath code2
addpath supp_material_RSCT/include

n_settings = 10;                  % # experimental settings

p_exp = [ 100 100 100 100 100 100 100 100 100 100];
d_exp = [  10  10  10  10  10  10  10  10  10  10];
L_exp = [   5   5   5   5   5   5  10  15  20  25];
n_exp = [  20  40  60  80 100  20  20  20  20  20];

n_algo = 5;                 % # algorithms : l1-min, OMP, Nuclear-norm min, Thresholding, NSN
n_trial = 10;               

T_exp = zeros(n_algo,n_settings);     % Computational time

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
    
    %% l1-min (SSC)
    i_algo = 1;
    tic;
    Z = admmLasso_mat_func(Y,false,20); Z_SSC = Z;
    T_exp(i_algo,i_exp) = T_exp(i_algo,i_exp) + toc;
    
    %% OMP (SSC-OMP)
    i_algo = 2; K = d-1;
    tic;
    Z = OMPSC(Y,K); Z_GFS = Z;
    T_exp(i_algo,i_exp) = T_exp(i_algo,i_exp) + toc;
      
    %% Nuclear norm minimization (LRR)
    i_algo = 3;
    tic;
    Z = solve_lrr(Y,1e5); Z_LRR = Z;
    T_exp(i_algo,i_exp) = T_exp(i_algo,i_exp) + toc;
    
    %% Thresholding (TSC)
    i_algo = 4; K = d-1;
    tic;
    Z = zeros(N,N);
    % the following scripts are copied from the TSC source code
    for i=1:N
        corvec = abs(Y'*Y(:,i));
        corvec(i) = 0; % so TSC will not select it
        [el,order] = sort(corvec, 'descend');
        Z(i, order(1:K) ) = exp(-2*acos(el(1:K))); % better than squared arcsin
    end
    Z = Z + Z';
    T_exp(i_algo,i_exp) = T_exp(i_algo,i_exp) + toc;

    %% NSN
    i_algo = 5; K = d-1;
    tic;
    [Z,~] = NSN(Y,K); Z_NSN = Z;
    T_exp(i_algo,i_exp) = T_exp(i_algo,i_exp) + toc;
    
    %%
    T_exp
end
end

close all;
n_algo = 5; algo = {'l1-minimization (SSC)', 'OMP (SSC-OMP)', 'Nuclear norm min. (LRR)', 'Thresholding (TSC)', 'NSN'};
clr = 'bgrcm'; mkr = 'osxd^'; 
figure;
set(gcf, 'Position', [100 100 1000 200]);
subplot(1,2,1);
for i=1:5
    plot(20:20:100,T_exp(i,1:5)/n_trial,[clr(i) mkr(i) '-']); hold on;
end
grid on; xlabel('Number of data points per subspace (n)'); ylabel('Time (sec)'); legend(algo);
subplot(1,2,2);
for i=1:5
    plot(5:5:25,T_exp(i,6:10)/n_trial,[clr(i) mkr(i) '-']); hold on;
end
grid on;xlabel('Number of subspaces (L)'); ylabel('Time (sec)'); legend(algo);

orient landscape;
set(gcf,'PaperPosition',[.25 3.0 10.5 2.5]);
print -dpdf test3


