function [Z,E] = solve_lrr23(X,lambda)
Q = orth(X');
A = X*Q;
[Z,E] = nlrr_dualS23a(X,A,lambda);
Z = Q*Z;