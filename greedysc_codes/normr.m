function Xnorm = normr(X)
    Xnorm = X./repmat(sqrt(sum(X.^2,2)),1,size(X,2));
end