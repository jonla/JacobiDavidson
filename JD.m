




function [lambda, e, res_hist, theta_approximations, count, V]=JD(A,guess)

%%% INITIALIZE:   Sets tolerance levels for convergence,
%%%              weight for the orthogonality of u to t,
%%%               normalizes the starting vector, calculates
%%%               initial eigenvalue approximation and declares
%%%               variables.
tol=10000;
tol2=0.00000001;
weight=1;
count=0;
iterations=length(guess);
I=eye(length(guess));
v1=guess/norm(guess);
V=[v1];
theta_init=v1'*A*v1;
orthogonality_V=[];
GUIDED=0;
NOT_GUIDED=0;
M=[theta_init];


theta_approximations = [];
res_hist = [];
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
    [eigvec, eigval] = eig(M);
    
    overlap = zeros(k,1);
    for j=1:k
        overlap(j) = abs(eigvec(1,j));
    end
    
    [~,indmax] = max(overlap);
    theta = eigval(indmax,indmax);
    s = eigvec(:,indmax)/norm(eigvec(:,indmax));
        
    
    theta_approximations = [theta_approximations theta];
    
    %%%% Calculates Ritz pairs and residual vector
    u=V*s;
    res=A*u-theta*u;
    
    res_hist = [res_hist norm(res)];
    
    %%% CORRECTION EQUATIONS
    if norm(res) > tol
        
        % guided correction               
        
        OCC=((I-u*u')*(A-theta_init*I)*(I-u*u'));
        %OCC_constrained=[OCC; weight*u'];
        %res_conc=[res; 0];
        t=OCC\-res;
        orthogonal_check=u'*t;   
               
    else if norm(res) > tol2
            
            % orthogonal completion correction
            
            OCC=((I-u*u')*(A-theta*I)*(I-u*u'));
            OCC_constrained=[OCC; weight*u'];
            res_conc=[-res; 0];
            tStart = tic;
            t=mldivide(OCC_constrained,res_conc);
            mld_time = toc(tStart)
            orthogonal_check=u'*t;
            
        else
            % Returns
            lambda=theta;
            e=u;
            break
        end        
    end
    
    %%% orthogonalize
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
    orthogonality_V=[orthogonality_V V(:,1)'*V(:,1+k)];
    
    
    %%% Constructs M
    for m=1:k
        M(m,k+1) = V(:,m)'*A*V(:,k+1);
        M(k+1,m) = V(:,k+1)'*A*V(:,m);
    end
    M(k+1,k+1) = V(:,k+1)'*A*V(:,k+1);
    
    
    
    
    
    count=count+1;
    end
        
end

