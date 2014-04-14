clc

k_max = 10;
tol = 1e-6;
dim = size(H,1);

which_state = 5;  % 297 for Hsmall
%guess = zeros(dim,1);
%guess(which_state) = (1 + 1i)/sqrt(2);
guess = ones(dim,1)/sqrt(dim);
lambda = guess'*H*guess;

[Q, R] = JD_harmritz(H,guess,lambda,k_max,tol);