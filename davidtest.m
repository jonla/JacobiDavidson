%% Symmetrisk matris storlek k
clc
k = 100;
A=[];
for i=1:k
    A(i,i)=i*10*(2*rand-1);
    %A(i,i)=100*complex(rand, 0.1*rand);
    for j=i+1:k
        A(j,i)=2*rand-1;
        %A(j,i)=complex(2*rand-1, 2*rand-1);
        A(i,j)=A(j,i);
    end
end

%%
clc
clf
data = load('singleparticle/SingleMatrixA_1.dat');
A = data(:,3)+1i*data(:,4);
clear data
k = sqrt(size(A, 1));
A = reshape(A, k, k);
imagesc(abs(A))
colormap(jet)

%%
clc
clf
tic
[v,e] = eig(A);
[~, i] = max(max(e));
emax = e(i,i);
vmax = v(:,i);
eigvalnr = 130;
e2=e(eigvalnr,eigvalnr);
v2=v(:,eigvalnr);
% v1 �r gissningen och �r n�stan = n�got egenv�rde
%v1 = v2 + 0.04*complex(2*rand(k,1)-1, 2*rand(k,1) -1);
%v1 = v2 + 0.08*(2*rand(k,1)-1);
v1 = zeros(k,1);
v1(5) = 1;
v1(81) = 1;
v1(82) = 1;
v1(83) = 1;
v1(84) = 1;

v1(86) = -1;
v1(87) = -1;
v1(88) = -1;
%v1 = ones(k,1);
v1 = v1/norm(v1);
accuracy = abs(v1'*v2)

subplot(2,1,2)
hold on
plot(diag(e),'*')
plot(e2, 'x', 'linewidth', 5)
%plot(sort(real(diag(e))), real(v2)/max(real(v2))*10^(-2))
plot(diag(e(125:160,125:160)),'*r')
toc
%%
%clf
clc
tic
m = 100;                  % Antal iterationer
V = v1;                   % Gissning
theta = v1'*A*v1;         % Approx till egenv�rde fr�n gissing
theta1 = theta;
subplot(2,1,2)
if imag(theta1) > 0
    plot(theta1,'xr')
else
    plot(0,theta1,'xr')
end
M = theta;
r = A*v1 - theta*v1;      % Residual
u = v1;
tmquad=[];
tgmres=[];
overlap=[];
orthcheck=[];
tKtildeinv= [];
tAtildeinv= [];
resnorm = [];

for i=1:m
    i
    % Davidson:
    % precond.: (diag(A)-Theta*I)
    %t = (diag(diag(A))-theta*eye(size(A)))\(-r);
    %t = t/norm(t);
    
    % Jacobi-Davidson correction EQ (t orth to u):
    %t1=toc;
    Atilde = (eye(size(A))-u*u')*(A-theta*eye(size(A)))*(eye(size(A))-u*u');
    %Ktilde = (eye(size(A))-u*u')*(diag(diag(A))-theta*eye(size(A)))*(eye(size(A))-u*u');
    LH = [Atilde; 100*u'];
    t = LH\[-r; 0];
    t = t/norm(t);
    
    % K=(diag(diag(A))-theta*eye(size(A)));
    % K=(A-theta*eye(size(A)));
    % Kinv=inv(K);
    % epsilon=(u'*Kinv*r)/(u'*Kinv*u);
    % t = epsilon*Kinv*u-Kinv*r;
    
    %t = t/norm(t);
    %t2=toc;
    %tmquad=[tmquad;t2-t1];              % time spent on min-quad and JD CEQ
    
    % Direct correctionequation
    %t1=toc;
    %tgm = dogmres(Atilde, -r, 5);
    %tgm = tgm/norm(tgm);
    %t2=toc;
    %tgmres=[tgmres;t2-t1]; % time spent on GMRES
    
    % GMRES with preconditioning
    t1=toc;
    K=diag(diag(A))-theta*eye(size(A));
    %K=A-theta*eye(size(A)) +0.4*(rand(k,k) - 0.5);
    %Ktilde = (eye(size(A))-u*u')*K*(eye(size(A))-u*u');
    %cond(inv(Ktilde)*Atilde)
    tgm = leftprecongmres(A, theta, K, r, u, 10);
    resnorm=[resnorm; norm(Atilde*tgm+r)];
    tgm = tgm/norm(tgm);
    t2=toc;
    tgmres=[tgmres;t2-t1]; % time spent on GMRES

    
    overlap=[overlap; t'*tgm];
    %orthcheck=[orthcheck; u'*tgm];
    t=tgm;
    
    % mod. Gram-Shmidt
    for j=1:i
        t = t - (V(:,j)'*t)*V(:,j);
    end
    vplus = t/norm(t);
    
    M = [M, (V'*A*vplus)
         (vplus'*A*V), vplus'*A*vplus];
    V = [V, vplus];
    [a,theta] = eig(M);
    %[~, p] = max(diag(abs(theta)));                % eigenval
    [~, p] = max(abs((a(1,:)))); 
    theta = theta(p,p);
    s = a(:,p);                                 % eigenvector to M
    u = V*s;                                    % Ritz vector (->eigenvec)
    r = A*u-theta*u;                          % Residual
    
    % Just for plotting
    subplot(2,1,1)
    hold on
    plot(i, log10(norm(r)),'*')
    ylabel('log10(norm(r))')
    grid on
    subplot(2,1,2)
    hold on
    plot(theta, 'or')
    ylabel('Eigenvalue')
    grid on
    
    % Check for convergance
    if norm(r) < 10^(-10)
        break;
    end
end
ttot=toc                        % Total time spent

subplot(2,1,2)
plot(theta, 'rx', 'linewidth', 2)
axis([0 0.4 -0.01 0.003])
% (theta, u) is best approx for eigenpair
% (theta1, v1) is first approx
% (e2, v2) is correct largest eigenpair (with eig(A))
time_ratio = sum(tgmres)/ttot
%%
clf
clc
subplot(2,1,1)
plot(abs(v2))
subplot(2,1,2)
plot(angle(v2))



%%
clf
hold on
K= @(omega) 1/(2-omega)*((1/omega)*(diag(diag(A))-tril(A,-1)))* ...
    (omega*inv(diag(diag(A))))*((1/omega)*(diag(diag(A))-triu(A,1)));
for i=1:1000
    Kn=norm(A-K(i*0.07+160));
    plot(i*0.07+160,Kn)
end

%%
tic
clc
clf
n = 25;
v1 = ones(k,1);
thety = [];
for i=1:k
    % v1 = rand([k,1]);
    [thetas, u] = davidson(A, v1, n);
    A = (eye(k)-u*u')*A*(eye(k)-u*u');
    if length(thetas)<n
        for j=1:n-length(thetas)
            thetas = [thetas, thetas(end)];
        end
    end
    thety = [thety; thetas];
    i
end
toc
plot(linspace(1,n+1,n+1),thety)
sort(thety(:,n))-sort(diag(e))
