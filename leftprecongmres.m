function [ x, x_history ] = leftprecongmres( A, theta, K, r, u, iter )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Parameters needed
%Kinv=inv(K);
uhat=K\u;
mu=u'*uhat;
rhat=K\r;
rtilde=rhat - (u'*rhat)/mu*uhat;

% Initiate GMRES
beta = norm(rtilde);
V = -rtilde/beta;
W = [];
H = [];
x_history = [];

for i=1:iter
    %W(:,i) = A*V(:,i);
    y = (A - theta*eye(size(A)))*V(:,i);
    yhat = K\y;
    W(:,i) = yhat - u'*yhat/mu*uhat;

    for j=1:i
        H(j,i) = W(:,i)'*V(:,j);
        W(:,i) = W(:,i) - H(j,i)*V(:,j);
        % W(:,i) = W(:,i) - (V(:,j)'*A*V(:,i))*V(:,j);
    end
    H(i+1,i) = norm(W(:,i)); % What if = 0?
    y = H\[beta; zeros(i,1)];
    x = V*y;
    x_history = [x_history x];
    V(:,i+1) = W(:,i)/H(i+1,i);
end
end

