function [ Q, R ] = JD_harmritz( ...
    A, guess, lambda, k_max, tol )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
dim = size(A,1);
m_max = 50;
m_min = 5;

%%% INITUAL PARAMETERS
t = guess;
k = 0;
m = 0;
Q = [];
R = [];
V = [];
VA = [];
W = [];
M = [];
MA = [];

while k <= k_max
    for i = 1:m
        t = t - (V(:,i)'*t)*V(:,i);
    end
    m = m + 1;
    V(:,m) = t/norm(t);
    VA(:,m) = A*V(:,m) - lambda*V(:,m);
    w = VA(:,m);
    
    for i = 1:k
        w = w - (Q(:,i)'*w)*Q(i,:);
    end

    
    for i = 1:(m-1)
        MA(i,m) = W(:,i)'*w;
        w = w - MA(i,m)*W(:,i);
    end
    MA(m,m) = norm(w);
    W(:,m) = w/MA(m,m);
    
    for i = 1:(m-1)
        M(i,m) = W(:,i)'*V(:,m);
        M(m,i) = W(:,m)'*V(:,i);
    end
    M(m,m) = W(:,m)'*V(:,m);
    
    % Generalized Schur:
    % SR and SL unitary,
    % TA and T upper triangular
    % SL*MA*SR = TA
    % SL*M*SR = T
    [TA, T, SR, SL] = qz(MA, M);
    SL = SL.';
    
    
    
    u = V*SR(:,1);
    uA = VA*SR(:,1);
    theta = T(1,1)'*TA(1,1);
    r = uA - theta*u;
    if size(Q,1) > 0
        atild = Q'*r;
        rtild = r - Q*atild;
    else
        atild = [];
        rtild = r;
    end
    
    while norm(rtild) < tol
        disp('Eigval found')
        R = [R atild
            zeros(size(atild))' theta + lambda];
        Q = [Q u];
        k = k + 1;
        if k == k_max return; end
        
        m = m - 1;
        for i = 1:m
            V(:,i) = V*SR(:,i+1);
            VA(:,i) = VA*SR(:,i+1);
            W(:,i) = W*SL(:,i+1);
            SR(:,i) = [zeros(i-1,1) 1 zeros(m-i+2)]'; % ??
            SL(:,i) =SR(:,i);
        end
        MA = TA(end-m+1:end,end-m+1:end);
        M = T(end-m+1:end,end-m+1:end);
        u = V(:,1);
        uA = VA(:,1);
        theta = M(1,1)'*MA(1,1);
        r = uA - theta*u;
        atild = Q*r;
        rtild = r - Q*atild;
    end
    
    % Restart procedures
    if m >= m_max
        disp('RESTART')
        for i = 2:m_min
            V(:,i) = V*SR(:,i);
            VA(:,i) = VA*SR(:,i);
            W(:,i) = W*SL(:,i);
        end
        MA = TA(1:m_min,1:m_min);
        M = T(1:m_min,1:m_min);
        V(:,1) = u;
        V = V(:,1:m_min);
        VA(:,1) = uA;
        VA = VA(:,1:m_min);
        W(:,1) = W*SL(:,1);
        W = W(:,1:m_min);
        m = m_min;
    end
    
    % CEQ with deflation
    THETA = theta + lambda;
    
    Qtild = [Q u];
    Kinv = sparse(1:dim,1:dim,1./(diag(A)-THETA));
    gmres_iter = 3;
    t = gmres_deflation(A,Kinv,Qtild,rtild,THETA,gmres_iter);
end
end

