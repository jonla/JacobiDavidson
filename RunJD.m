
%%% Full tester for Jacobi Davidson with GMRES solver for the CEQ %%%

%%%===============================================================%%%

%  
% %%% ==== Loads single particle matrices ==== %%%
% importer = load('/net/data2/riklund/ThesisTwoparticle1/SingleMatrixA_1.dat');
% A1 = importer(:,3)+1i*importer(:,4);
% N = sqrt(length(importer(:,1)));
% clear importer
% A1 = reshape(A1,N,N);
% 
% importer = load('/net/data2/riklund/ThesisTwoparticle1/SingleMatrixA_2.dat');
% A2 = importer(:,3)+1i*importer(:,4);
% N = sqrt(length(importer(:,1)));
% clear importer
% A2 = reshape(A2,N,N);
% 
% disp('Single particle matrices loaded')
% 
% %%% Forms matrix D
% [~, D_1] = eig(A1);
% [~, D_2] = eig(A2);
% 
% D = sparse([]);
% for i = 1 : N
%     D((i-1)*N+1:i*N, (i-1)*N+1:i*N) = D_2 + eye(size(D_2))*D_1(i,i);
% end
% disp('D constructed')
% 
% % Load A and construct H = A + D
% data=load('/net/data2/riklund/ThesisTwoparticle1/matrixA.dat');
% disp('A loaded')
% data(:,1:2) = [];
% H = data(:,1)+1i*data(:,2);
% clear data
% dim = sqrt(size(H,1));
% H = reshape(H,dim,dim);
% disp('Matrix A created')
% H = H + D;
% disp('Full matrix H created')
% clear D


%%% ===== Sets initial guess ===== %%%
%    Resonance states for sp matrices are 120 and 129. A suitable
%    guess could be 160*129+121 = 20761
dim=size(H,1);
%which_state=1300; 
guess=zeros(dim,1);
guess(which_state,1) = 1+1i;
init_lambda=guess'*H*guess;
disp('Guess vector and initial theta calculated')


disp('Running JD...')
%%%  ===== Run JD with gmres ===== %%%
clc
diary('Timelogs.txt')
format compact
tStart_gm = tic ;  
[lambda, state, reshist, theta_hist, count] = davidson(H,guess,1e-10,'T');
JD_Time = toc(tStart_gm)

disp('JD complete!')

disp('Found eigenvalue:')
lambda
diary end
