function Z = Run_SpNM(method,Xo,lambda)

switch method
     case 'spdualNulrr'
        [Z,E] = solve_lrr(Xo,lambda);
        
    case 'spdual23lrr'
        [Z,E] = solve_lrr23(Xo,lambda);
        
    case 'spdual12lrr'
        [Z,E] = solve_lrr12(Xo,lambda);   
         
end

% 
% % ����Ĵ����MS�н��������� ��FC������������ 
% % for i = 1 : size(L,2)  
% %    L(:,i) = L(:,i) / max(abs(L(:,i))) ;    
% % end   % normalized coefficient matrix

%Z = processC(Z,0.9);

% refining Z
%[U,S,V] = svd(Z);
%S = diag(S);
%r = min(4*K+1,sum(S>1e-3*S(1)));
%S = S(1:r);
%U = U(:,1:r)*diag(sqrt(S));
%U = normr(U);
%Z = U*U';Z=abs(Z);
%L = Z.^4;

% % spectral clustering
% L = (L + L')/2;
% D = diag(1./sqrt(sum(L,2)));
% L = D*L*D;
% CKSym = L;

% nCluster = length(unique(gnd));
% idx = clu_ncut(L,nCluster);  % turnable 
% acc = compacc(idx,gnd);  
%nCluster = max(s);
%nCluster = length(unique(gnd));
%idx = clu_ncut(L,nCluster);  % turnable
%
%missrate = 1-compacc(idx',gnd);

end