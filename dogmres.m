function [ x ] = dogmres( A, b, theta, K, u, iterations )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
create_rtilde_time = 0;
create_y_time = 0;
build_H_time = 0;
solve_y_time = 0;
build_wk_time = 0;
I = sparse(eye(size(A)));
mu = K\u;
t1 = tic;
rhat = K\b;
rtilde = rhat-u'*rhat*mu/(u'*mu);
create_rtilde_time = create_rtilde_time + toc(t1)
beta = norm(rtilde);
V = rtilde/beta;
W = [];
H = [];
x_history = [];

for i=1:iterations
    t2 = tic;
    y = (A-theta*I)*V(:,i);
    create_y_time =create_y_time + toc(t2)
    t5 = tic;
    yhat = K\y;
    W(:,i) = yhat-u'*yhat/(u'*mu)*mu;
    build_wk_time = build_wk_time + toc(t5)
    t3 = tic;
    for j=1:i
        H(j,i) = W(:,i)'*V(:,j);
        W(:,i) = W(:,i) - H(j,i)*V(:,j);
        
        
    end
    wk_norm = norm(W(:,i));
    if wk_norm < 10^-10
        
        x = V*y;
        break
    end
    H(i+1,i) = wk_norm;
    build_H_time = build_H_time + toc(t3)
    
    t4 = tic;
    y = H\[beta; zeros(i,1)];
    solve_y_time = solve_y_time + toc(t4)
    %[U,Sigma,Q] = svd(H);
    %y = 0;
    % for k = 1:i-1
    %     y = y+Q(:,k)*U(:,k)'*[beta; zeros(i,1)]/Sigma(k,k);
    % end
    
    
    
    
    V(:,i+1) = W(:,i)/H(i+1,i);
end
x = V(:,1:end-1)*y;

