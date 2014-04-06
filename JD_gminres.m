




function [lambda, e, res_hist, theta_approximations, count]=JD_gminres(A,guess,GmresIterations,tol2,Time)

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
    disp('New iteration beginning')
    
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
    disp('M solved for eigenpair')
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
    disp('Solving CEQ...')
    if norm(res) > tol
        
        % guided correction
        disp('Guided with initial theta')
        K = sparse(diag(diag(A))-theta*I);
        t = dogmres(A,-res,theta,K,u,GmresIterations);
        
        
        
        
        
    else if r_norm > tol2
            
            disp('Guided with updated theta')
            
            
            tStart = tic;
            K = sparse(diag(diag(A))-theta*I);
            t = dogmres(A,-res,theta,K,u,GmresIterations);
            %t = gmres(A,-res,[],0.1,5,Ktilde);
            gmres_time = gmres_time + toc(tStart);
            GMRES_TIME = [k gmres_time];
            
            
            
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
        disp('Done!')
        %%% orthogonalize
        disp('Orthogonalzing...')
        if mod(k,15) == 0
            disp('Full reorthogonalization...')
            for i=1:k-1
                for j=1:i
                    V(:,i+1) = V(:,i+1)-V(:,j)'*V(:,i+1)*V(:,j)
                end
            end
        else
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
        end
        disp('Done!')
        %%% Expands search space
        v=t/norm(t);
        V=[V v];
        gs_time = gs_time + toc(t1);
        GS_TIME = [k gs_time];
        
        
        
        %%% Constructs M
        % This part could be improved by preallocating the size of M
        t2 = tic;
        disp('Building M...')
        for m=1:k
            M(m,k+1) = V(:,m)'*A*V(:,k+1);
            M(k+1,m) = V(:,k+1)'*A*V(:,m);
        end
        M(k+1,k+1) = V(:,k+1)'*A*V(:,k+1);
        M_time = M_time + toc(t2);
        M_TIME = [k M_time];
        disp('Done!')
        
        
        count=count+1;
    end
    
end

