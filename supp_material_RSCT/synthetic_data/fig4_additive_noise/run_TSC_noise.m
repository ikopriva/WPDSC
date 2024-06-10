
%%%%%%%%%%%%%%%%%%
% performance of TSC under additive noise 
%%%%%%%%%%%%%%%%%%

function run_TSC_noise(L,m,nit,dvec,nl,sig,filename)

addpath '../../include'

ce = zeros(length(dvec),length(nl));

for j=1:length(dvec)  % vary over d
	j    
    for k=1:length(nl)% vary over nl
		nl(k)        
        N = L*nl(k);  % number of points overall

        for i=1:nit
            
			X = [];
            Xlabels = [];
            for l=1:L
                U = orth(randn(m,dvec(j)));
                A = normc(randn(dvec(j),nl(k))); 
                X = [X, U*A];  
                Xlabels = [Xlabels, ones(1,nl(k))*l];
            end            
            
           	X = X + sig/sqrt(m)*randn(m,N); 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            q = 2*max(3,ceil(nl(k)/20));
            [labels,A,Lest] = TSC(X,q);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % compute error 
            ce(j,k) = ce(j,k) + clustering_error(labels',Xlabels)/nit;
            
		end
        
    end
end

em2tikzf(ce,dvec,nl,filename);

end
