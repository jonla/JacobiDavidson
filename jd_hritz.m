




function [lambda, e, res_hist, theta_approximations, count]=JD_gminres(A,guess,GmresIterations,tol2,Time)

%%% INITIALIZE:   Sets tolerance levels for convergence,
%%%              weight for the orthogonality of u to t,
%%%               normalizes the starting vector, calculates
%%%               initial eigenvalue approximation and declares
%%%               variables.
tol=10000;
maxit = 5;
%tol2=10^-6;
dim=length(guess);
I=sparse(eye(dim));
v1=guess/norm(guess);
V=[v1];
W=[A*v1/norm(A*v1)];
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


for m = 1:maxit
    
    %%%% INNER LOOP %%%%
    for k=1:10
        disp('New inner iteration beginning')
        
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
        u=V*s/norm(V*s);
        uhat=A*u;
        res=uhat-theta*u;
        disp('Residual norm:')
        r_norm = norm(res)
        res_hist = [res_hist r_norm];
        Ritz_time = Ritz_time + toc(t4);
        RITZ_TIME = [k Ritz_time];
        
        %%% CORRECTION EQUATIONS %%%
        disp('Solving CEQ...')
        if norm(res) > tol
                       
            disp('Guided with initial theta')
            K = sparse(diag(diag(A))-theta*I);
            t = dogmres(A,-res,theta,K,u,GmresIterations);
                       
        else if r_norm > tol2
                
                disp('Guided with updated theta')
                tStart = tic;
                Kinv = sparse(1:dim,1:dim,1./(diag(A)-theta));
                t = harm_gmres(A,-res,theta,Kinv,u,uhat,GmresIterations);
                %t=Kinv*res;
                gmres_time = gmres_time + toc(tStart);
                GMRES_TIME = [k gmres_time];
                
            else
                %%% Returns %%%
                lambda=theta;
                e=u;
                if strcmp(Time,'T')
                    disp('Timings:')
                    M_EIG_ = M_EIG_TIME
                    RITZ = RITZ_TIME
                    GMRES = GMRES_TIME
                    GS = GS_TIME
                    M_BUILD = M_TIME
                    break
                end
            end
        end
        disp('Done!')
        
        
        %%% orthogonalize V %%%
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
                t=t-V(:,i)*W(:,i)'*t;
                if norm(t)/norm(t_prim) < 0.250
                    for j=1:k
                        t=t-V(:,j)*W(:,j)'*t;
                    end
                end
            end
        end
        disp('Done!')
        
        
        %%% Expands search space V
        v=t/norm(W(:,i));
        V=[V v];
        gs_time = gs_time + toc(t1);
        GS_TIME = [k gs_time];
        
        %%% Orthogonalize and expand W
        w = A*v;
        for i=1:k
            w_prim = w;
            w=w-W(:,i)'*w*W(:,i);
            if norm(w)/norm(w_prim) < 0.250
                for j=1:k
                    w=w-W(:,i)*W(:,i)'*w;
                end
            end
        end
        w=w/norm(w);
        W = [W w];
        
        %%% Constructs M on the orthogonal space of AV
        % This part could be improved by preallocating the size of M
        t2 = tic;
        disp('Building M...')
       
        
        for m=1:k
            M(m,k+1) = W(:,m)'*V(:,k+1);
            M(k+1,m) = W(:,k+1)'*V(:,m);
        end
        M(k+1,k+1) = W(:,k+1)'*V(:,k+1);
        M_time = M_time + toc(t2);
        M_TIME = [k M_time];
        M = inv(M);
        
        
        count=count+1
    end
    %%% Resetting spaces and variables
    disp('RESTARTING...')
    M = [theta];
    V = [u];
end

end

