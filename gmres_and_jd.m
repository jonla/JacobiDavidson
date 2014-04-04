




function [lambda, e, res_hist, theta_approximations, count]=JDgminres(A,guess,GmresIterations,tol2,Time)

%%% INITIALIZE:   Sets tolerance levels for convergence,
%%%              weight for the orthogonality of u to t,
%%%               normalizes the starting vector, calculates
%%%               initial eigenvalue approximation and declares
%%%               variables.
tol=10000;
%tol2=10^-6;
iterations=length(guess);
I=sparse(eye(length(guess)));
v1=guess/norm(guess);
V=[v1];
theta_init=v1'*A*v1;



M=[theta_init];
theta_approximations = [];
res_hist = [];
count=0;


%%% Timings
gmres_time = 0;
gs_time = 0;
M_time = 0;
Meig_time = 0;
Ritz_time = 0;



%%%% OUTER LOOP
for k=1:iterations
        
    
    %%% This part finds eigenpair closest to guess
%     [eigvec,eigval]=eig(M);      
%     diff=diag(eigval)-theta_init;
%     [~, ind] = min(abs(diff));
%     theta=eigval(ind,ind);
%     s=eigvec(:,ind)/norm(eigvec(:,ind))
    
    
    %%% Finds largest eigenpair of M
%     [eigvec,eigval] = eig(M);
%     [~,indmax] = max(abs(diag(eigval)));
%     theta = eigval(indmax,indmax);
%     s = eigvec(:,indmax)/norm(eigvec(:,indmax));
    
    %%% Finds eigenpar with largest overlapp to v0
    t3 = tic;
    [eigvec, eigval] = eig(M);
    
    overlap = zeros(k,1);
    for j=1:k
        overlap(j) = abs(eigvec(1,j));
    end
    
    
    Meig_time = Meig_time + toc(t3);
    M_EIG_TIME = [k Meig_time];
    
    t4 = tic;
    [~,indmax] = max(overlap);
    theta = eigval(indmax,indmax);
    s = eigvec(:,indmax)/norm(eigvec(:,indmax));
        
    
    theta_approximations = [theta_approximations theta];
    
    %%%% Calculates Ritz pairs and residual vector
    u=V*s;
    res=A*u-theta*u;
    disp('Residual norm:')
    r_norm = norm(res)
    res_hist = [res_hist r_norm];
    Ritz_time = Ritz_time + toc(t4);
    RITZ_TIME = [k Ritz_time];
    
    %%% CORRECTION EQUATIONS
    if norm(res) > tol
        
        % guided correction               
        
        OCC=((I-u*u')*(A-theta_init*I)*(I-u*u'));
        
        
        
        
               
    else if r_norm > tol2
            % orthogonal completion correction
            
            %%% MLDIVIDE                    
%             OCC=((I-u*u')*(A-theta*I)*(I-u*u'));
%             OCC_constrained=[OCC; weight*u'];
%             res_conc=[-res; 0];
%             t_mld=mldivide(OCC_constrained,res_conc);
%             mldNorm = norm(t_mld);
            

            %%% GMRES
                
                %t1 = tic;
                %Ktilde = ((I-u*u')*K*(I-u*u'));
                %tilde_time = tilde_time + toc(t1)
                tStart = tic;
                K = sparse(diag(diag(A))-theta*I);
                t = dogmres(A,-res,theta,K,u,GmresIterations);
                %t = gmres(A,-res,[],0.1,5,Ktilde);
                gmres_time = gmres_time + toc(tStart);
                GMRES_TIME = [k gmres_time];
                
            %%% Comparisons
                %dogNorm = norm(t);
                %orthogonal_check=u'*t;
                %cond_Atilde = cond(OCC)
                %cond_prec = cond(inv(Ktilde)*OCC)
                %t=t/norm(t);
                %t_mld = t/norm(t_mld);
                %t_overlap = [t_overlap; abs(t_mld'*t)];
            
        else
            % Returns
            lambda=theta;
            e=u;
            if strcomp(Time,'T')
                disp('Timings:')
                M_EIG_ = M_EIG_TIME
                RITZ = RITZ_TIME
                GMRES = GMRES_TIME
                GS = GS_TIME
                M_BUILD = M_TIME
            break
        end        
    end
    
    %%% orthogonalize
    t1 = tic;
    for i=1:k 
        t_prim=t;
        t=t-V(:,i)'*t*V(:,i);
        if norm(t)/norm(t_prim) < 0.250
            for j=1:k
                t=t-V(:,j)'*t*V(:,j);
            end
        end
    end
    
    %%% Expands search space
    v=t/norm(t);
    V=[V v];
    gs_time = gs_time + toc(t1);
    GS_TIME = [k gs_time];
    
    
    
    %%% Constructs M
    % This part could be improved by preallocating the size of M
    t2 = tic;
    for m=1:k
        M(m,k+1) = V(:,m)'*A*V(:,k+1);
        M(k+1,m) = V(:,k+1)'*A*V(:,m);
    end
    M(k+1,k+1) = V(:,k+1)'*A*V(:,k+1);
    M_time = M_time + toc(t2);
    M_TIME = [k M_time];
    
    
    
    count=count+1;
    end
        
end

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

