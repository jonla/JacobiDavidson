function [lambda, e, res_hist, theta_hist, count]=davidson(A,guess,tol2,Time)

%%% INITIALIZE:   Sets initial parameters needed
%%%               for first outer loop.

dim = size(A,1);
iterations = 300;  % Never do more then 300 iterations
%I = sparse(eye(length(guess)));
v1 = guess/norm(guess);
V = v1;
theta_init = v1'*A*v1;

M = theta_init;
theta_hist = [];
res_hist = [];
count = 0;

%%% Timings
ceq_time = 0;
gs_time = 0;
M_time = 0;
Meig_time = 0;
Ritz_time = 0;

%%%% OUTER LOOP
for k = 1:iterations
    %%% Finds eigenpar with largest overlapp to v0
    t3 = tic;
    [eigvec, eigval] = eig(M);
    [~, indmax] = max(abs((eigvec(1,:))));
    %[~, indmax] = min(abs(diag(eigval)-theta_init));
    Meig_time = Meig_time + toc(t3);
    
    t4 = tic;
    theta = eigval(indmax,indmax);
    s = eigvec(:,indmax)/norm(eigvec(:,indmax));
    %disp(['Eigval. it ', num2str(k), ': ', num2str(theta)])
    
    theta_hist = [theta_hist theta];
    
    %%%% Calculates Ritz pairs and residual vector
    u = V*s;
    res = A*u - theta*u;
    r_norm = norm(res);
    disp(['Residual norm: ', num2str(r_norm)])
    res_hist = [res_hist r_norm];
    Ritz_time = Ritz_time + toc(t4);
    RITZ_TIME = [k Ritz_time];

    if r_norm > tol2
        %%% Davidsons CEQ
        tStart = tic;
        Kinv = sparse(1:dim,1:dim,1./(diag(A)-theta));
        %t = dogmres(A,-res,theta,Kinv,u,GmresIterations);
        t = Kinv*res;

        ceq_time = ceq_time + toc(tStart);
    else
        % Convergence -> Returns
        lambda = theta;
        e = u;
        if strcmp(Time,'T')
            disp('Timings:')
            RITZ = RITZ_TIME
            CEQ = [k, ceq_time]
            GS = GS_TIME
            M_BUILD = M_TIME
        end
        break
    end
    
    %%% orthogonalize
    t1 = tic;
    for i = 1:k 
        t_prim = t;
        t = t - V(:,i)'*t*V(:,i);
        if norm(t)/norm(t_prim) < 0.250
            for j = 1:k
                t = t - V(:,j)'*t*V(:,j);
            end
        end
    end
    
    %%% Expands search space
    v = t/norm(t);
    V = [V v];
    gs_time = gs_time + toc(t1);
    GS_TIME = [k gs_time];
    
    
    
    %%% Constructs M
    % This part could be improved by preallocating the size of M
    t2 = tic;
    W1 = A*(V(:,k+1));
    W2 = V(:,k+1)'*A;
    for m = 1:k
        M(m,k+1) = V(:,m)'*W1;
        M(k+1,m) = W2*V(:,m);
    end
    M(k+1,k+1) = V(:,k+1)'*W1;
    M_time = M_time + toc(t2);
    M_TIME = [k M_time];
    
    
    
    count = count + 1;
end
        
end

