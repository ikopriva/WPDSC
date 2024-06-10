clear all; close all; clc;

% Before running this code, SSC and LRR codes from the authors' websites
% should be in the subfolders with the following titles
addpath SSC_ADMM_v1.1
addpath code2
addpath supp_material_RSCT/include
addpath SCC

load Hopkins155_titles.mat; i_data = 0;
for i=1:length(Hopkins155_titles)
    fprintf('%d/%d %s \n',i,length(Hopkins155_titles),Hopkins155_titles{i});
    eval(['load Hopkins155/' Hopkins155_titles{i} '/' Hopkins155_titles{i} '_truth.mat' ]);
    
    L = max(s);

    if (L == 3)
        
        i_data = i_data + 1;
        
        N = size(x,2);
        F = size(x,3);
        p = 2*F;
        Y = reshape(permute(x(1:2,:,:),[1 3 2]),p,N);
        
        A0 = s; [~,I] = sort(A0,1);

%% K-means
        fprintf('Running K-means..\n'); i_algo = 1;
        
        tic;
        [A,~] = kmeans(Y',L,'emptyaction','singleton','replicates',10,'display','off');
        ET(i_algo,i_data)  = toc;
        CE(i_algo,i_data)  = computeCE(A,A0);
        
        
%% K-flats
        fprintf('Running K-flats..\n'); i_algo = 2;
        
        CE(i_algo,i_data)  = 1;
        for i_repl = 1:10
            tic;
            A = Kflats(Y,3,L);
            ET(i_algo,i_data)  = toc * 10;
            CE(i_algo,i_data)  = min(computeCE(A,A0), CE(i_algo,i_data));
        end

%% SSC
        fprintf('Running SSC..\n'); i_algo = 3;
        r = 0; affine = true; outlier = false; rho = 0.7; alpha = 800;
        tic
        [~,Z,A] = SSC(Y,r,affine,alpha,outlier,rho,s);
        
        ET(i_algo,i_data)  = toc;
        CE(i_algo,i_data)  = computeCE(A,A0);
       
%% LRR
        fprintf('Running LRR..\n'); i_algo = 4;
        % The following scripts are copied from the LRR source code.
        tic
        lambda = 4;
        %run lrr
        Z = solve_lrr(Y,lambda);
        %post processing
        [U,S,V] = svd(Z,'econ');
        S = diag(S);
        r = sum(S>1e-4*S(1));
        U = U(:,1:r);S = S(1:r);
        U = U*diag(sqrt(S));
        %U = normr(U);
        U = U./repmat(sqrt(sum(U.^2,2)),1,size(U,2));
        LL = (U*U').^4;
        % spectral clustering
        D = diag(1./sqrt(sum(LL,2)));
        LL = D*LL*D;
        [U,S,V] = svd(LL);
        V = U(:,1:L);
        V = D*V;
        A = kmeans(V,L,'emptyaction','singleton','replicates',10,'display','off');
        
        ET(i_algo,i_data)  = toc;
        CE(i_algo,i_data)  = computeCE(A,A0);

%% SCC
        fprintf('Running SCC..\n'); i_algo = 5;
        tic
        [A,~] = scc(Y',3,L);
        ET(i_algo,i_data)  = toc;
        CE(i_algo,i_data)  = computeCE(A,A0);        

%% GFS
        fprintf('Running GFS..\n'); i_algo = 6;
        tic
        Z = OMPSC(Y,8);
        A = SpectralClustering(abs(Z)+abs(Z)',L);
        ET(i_algo,i_data)  = toc;
        CE(i_algo,i_data)  = computeCE(A,A0);

%% TSC
        fprintf('Running TSC..\n'); i_algo = 7;
        tic
        [A,Z] = TSC(Y,10,L);
        ET(i_algo,i_data)  = toc;
        CE(i_algo,i_data)  = computeCE(A,A0);
         
%% NSN+Spectral
        fprintf('Running NSN+Spectral..\n'); i_algo = 8;
        tic
        Y = Y./repmat(sqrt(sum(Y.^2,1)),p,1);
        Z = NSN(Y,5,5,1e-4);
        A = SpectralClustering(Z+Z',L);
            
        ET(i_algo,i_data)  = toc;
        CE(i_algo,i_data)  = computeCE(A,A0);
        
        CE(:,i_data)
%%               
    end
    
end
            mean(CE')
            median(CE')
            max(CE')
            mean(ET')    



