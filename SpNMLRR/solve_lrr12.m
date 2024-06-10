function [Z,E] = solve_lrr12(X,lambda)
Q = orth(X');
A = X*Q;
[Z,E] = nlrr_dualS12a(X,A,lambda);
Z = Q*Z;