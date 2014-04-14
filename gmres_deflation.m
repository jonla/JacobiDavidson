function [ x ] = gmres_deflation( A,Kinv,Qtild,r,theta,gmres_iter )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

%%% INITALIZE
%dim = length(r);
Qhat = Kinv*Qtild;
M = Qtild'*Qhat;
[L, U] = lu(M);

rhat = Kinv*r;
gamma = Qtild'*rhat;
betha = L\gamma;
alpha = U\betha;
rtild = rhat - Qhat*alpha;

%%% GMRES
beta = norm(rtild);
V = -rtild/beta;
W = [];
H = [];

for i=1:gmres_iter
    y = A*V(:,i) - theta*V(:,i);
    yhat = Kinv*y;
    
    gamma = Qtild'*yhat;
    betha = L\gamma;
    alpha = U\betha;
    W(:,i) = yhat - Qhat*alpha;
    
    for j=1:i
        H(j,i) = W(:,i)'*V(:,j);
        W(:,i) = W(:,i) - H(j,i)*V(:,j);
    end
    
    H(i+1,i) = norm(W(:,i)); % What if = 0?
    y = H\[beta; zeros(i,1)];
    x = V*y;
    V(:,i+1) = W(:,i)/H(i+1,i);
end
x = x/norm(x);
end

