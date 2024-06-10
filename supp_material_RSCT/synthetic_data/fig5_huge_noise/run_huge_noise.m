
%%%%%%%%%%%%%%%%%%
% performance of TSC under huge additive noise 
%%%%%%%%%%%%%%%%%%


% original parameters
L = 5;
m = 400;
nit = 10;
nl = 10:10:100;
d = 5;
sigvec = 0:0.5:4;

% parameters to produce a figure fast
L = 2;
m = 400;
nit = 1;
nl = 10:40:100;
d = 5;
sigvec = 0:2:4;

addpath '../../include'

ce = zeros(length(sigvec),length(nl));
fde = zeros(length(sigvec),length(nl));
Le = zeros(length(sigvec),length(nl));
Ls = zeros(length(sigvec),length(nl));

for j=1:length(sigvec)  % vary over sigma
    
    for k=1:length(nl)% vary over nl
        
        N = L*nl(k);  % number of points overall

        for i=1:nit
            X = [];
            Xlabels = [];
            for l=1:L
                U = orth(randn(m,d));
                A = normc(randn(d,nl(k))); 
                X = [X, U*A];  
                Xlabels = [Xlabels, ones(1,nl(k))*l];
            end            
            
           	X = X + sigvec(j)/sqrt(m)*randn(m,N); 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            q = 4*max(3,ceil(nl(k)/20));
            [labels,A,Lest] = TSC(X,q);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % compute error 
            ce(j,k) = ce(j,k) + clustering_error(labels',Xlabels)/nit;
           	fde(j,k) = fde(j,k) + feature_detection_error(A,L)/nit; 
		
	       	if(L ~= Lest)
                Le(j,k) = Le(j,k) + 1/nit;
                if L < Lest
                    Ls(j,k) = Ls(j,k) + 1/nit;
                else
                    Ls(j,k) = Ls(j,k) - 1/nit;
                end
            end	
		
		
		end
        
    end
end


em2tikzf(Le,sigvec,nl,'./LE_huge_noise.dat');
em2tikzf(Ls,sigvec,nl,'./LS_huge_noise.dat');
em2tikzf(fde,sigvec,nl,'./FDE_huge_noise.dat');
em2tikzf(ce,sigvec,nl,'./CE_huge_noise.dat');

