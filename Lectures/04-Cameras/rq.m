function [R,Q]= rq(A)
% full RQ factorization
% returns [R 0] and Q such that [R 0] Q = A
% where R is upper triangular and Q is unitary (square)
% and A is m by n, n >= m  (A is short and fat)
% Equivalent to   A' = Q' [R']  (A' is tall and skinny)
%                         [0 ]
% or    
%                 A'P = Q'[P 0] [P R' P]
%                         [0 I] [  0   ]
% where P is the reverse identity matrix of order m (small dim) 
% This is an ordinary QR factorization, say Q2 R2.
% Written by Michael Overton, overton@cs.nyu.edu (last updated Nov 2003)

m = size(A,1);
n = size(A,2);
if n < m
   error('RQ requires m <= n')
end
P = fliplr(eye(m));
AtP = A'*P;
[Q2,R2] = qr(AtP);
bigperm = [P zeros(m,n-m); zeros(n-m,m) eye(n-m)];
Q = (Q2*bigperm)';
R = (bigperm*R2*P)';
