clc
clear all

%% Builds test matrices

N=1000;
I = eye(N);
Breal=complex(-1+2*rand(N,N));
Bcomplex = i*(-1+2*rand(N,N));
A_org = Breal+Bcomplex;
A_hermitian=(A_org+A_org')/2;
A_complex = (A_org+transpose(A_org))/2;
Areal = (Breal+Breal')/2;
Asparse = zeros(N,N);
    D = 50; %%% Number of non zero elements on a row
    for i=1:N
        Asparse(i, i) = 2 *rand(1,1) - 1;
        for j=1:i-1
            if abs(i - j) < D
                Asparse(i, j:j + 1) = (2 * rand(1,1)-1+ 2 *rand(1,1) * 1j + 1j);
            end
                Asparse(j, i) = conj(Asparse(i, j));
        end
    end


%% Choose matrix type
A = Asparse;



%% Constructs initial vector
which_state=100;
tStart = tic;
[E,V]=eig(A);
tElapsed = toc(tStart)
eigs=diag(V);
sought_state=E(:,which_state);
sought_lambda=eigs(which_state);
error=0.015*ones(N,1);
guess=(sought_state+error)/norm(error+sought_state);
guess_accuracy=abs(sought_state'*guess);

init_lambda=guess'*A*guess;



cond_A = cond(A)
precond_cond = cond(inv(diag(diag(A))-init_lambda*I)*(A-init_lambda))



%% Run JD with mldivide. 
clc

tStart_lstsq = tic;
[lambda_guided, e_guided, res, prev_theta, count, V] = JD(A,guess);
tElapsed_FULL = toc(tStart_lstsq);
%% Run JD with gmres
clc
  tStart_gm = tic ;  
  [lambda_guidedgm, e_guidedgm, resgm, prev_thetagm, countgm, Vortgm] = JD_gminres(A,guess);
  tElapsed_GMRES = toc(tStart_gm);

%% Real matrix plot
clf
subplot(3,1,1)
plot(count+1,lambda_guided,'o')
hold on
plot(count+1,init_lambda,'og')
hold on
plot(count+1,sought_lambda,'*r')
plot(prev_theta,'k')
plot(count+1,eigs,'.k')
hold on
subplot(3,1,2)
plot(e_guided,'r')
hold on
plot(E(:,which_state))
subplot(3,1,3)
plot(log10(res))

%% Complex matrix plot
clf
clc
its = linspace(0,N,N);
subplot(2,1,1)
plot(real(prev_theta),imag(prev_theta),'r')
hold on
plot(real(eigs),imag(eigs),'.k')
hold on
plot(real(lambda_guided),imag(lambda_guided),'o')
hold on
plot(sought_lambda,'g*')
phase = exp(-1i*0.129);
e_phased = (e_guided'*sought_state)*e_guided;
subplot(2,1,2)
plot(log10(res))
grid on
shg
%% 3D-plots
clf
plot3(real(sought_state),imag(sought_state),its)

hold on
plot3(real(e_phased),imag(e_phased),its,'r')


shg