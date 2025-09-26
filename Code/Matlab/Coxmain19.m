%% Format
clc; clear; close all; format short;  % format long;
%% Global variables
global L1  X_train  ytime_train  R_matrix  ystatus_train  lambda_opt  alpha_opt  m  p1

%% set pathway
path = '/home/lly/R/CNetCox/Data/';
% path = '/Users/lilingyu/E/PhD/R/CNetCox/Data/';
% path = 'D:\E\博士\R_程序\CNetCox\Data\';
[path,'TCGA_NEW/']

%% Adjacency matrix
adj = importdata([path,'TCGA_NEW/adjmatrix_comp_UNG.csv']);
% adj = importdata('adjmatrix_comp_UNG.csv');
A_adj = sparse(adj.data);
%% Set folds in cv
nfold = 5;

%% Cut point
vector_hat = importdata([path,'TCGA_NEW/cut_vector_UNgene.txt']); 
% vector_hat = importdata('cut_vector_UNgene.txt'); 
delta_hat = vector_hat.data;

%% Set options
options = optimoptions('fmincon','Algorithm','interior-point');
options.MaxFunEvals = 1e6;

%% Inlut data
i = '19';
txt = importdata([path,'Data_train/',i,'.txt']);
% txt = importdata([path,'Data_train/',int2str(i)','.txt']); 
% txt = importdata(['Data_train/23.txt']);  
train_data = txt.data;
[m,q] = size(train_data);
% Cox need to begin at third col
X_train = train_data(1:m,3:q);      
ytime_train = train_data(1:m,1);     
ystatus_train = train_data(1:m,2);    
p = q-2;
p1= q-1;

%% laplace Matix
[L, L1] = Laplacian_Matrix(p,A_adj);
% L = eye(p,p);

%% Risk Matrix
R_matrix = RiskMatrix(ytime_train);
% spy(R_matrix)

%% candidate lambdas in simulations
% e=(log(100)-log(.0001))/99;
% lambda=exp(log(.0001):e:log(100)); 
%% set candidate values of alpha 
% alpha1=0.2:.2:.6;  
%% provide the value of alpha
alpha = 0.5;  

%% get the least upper bound, lambda_max
lammax = getLambMaxCox(X_train, ystatus_train, alpha); 
e = (log(lammax)-log(1))/19;
lambda = exp(log(1):e:log(lammax)); 

%% get the optimal lambda and alpha through cross validation 
[lambda_opt, alpha_opt, r] = cvCox( X_train, ytime_train, ystatus_train, R_matrix, L1, alpha, lambda, nfold );

%% Optimal parameters
alpha_opt 
lambda_opt 

%% Constraints
theta_0 = zeros(p1,1);    % Initial value is 0
u_0 = (1e-4)*ones(p1,1);  % The initial value is 1

a11 = eye(p1,p1);
a22 = -1*eye(p1,p1);
A1 = [a11;a22];
A2 = [a22;a22];
A = [A1,A2];
b = zeros(2*p1,1);

% A = [];
% b = [];
Aeq = []; 
beq = [];
% vlb = [];        
% vub = [];  

vlb1 = zeros(p1,1);        
vub1 = zeros(p1,1);  
delta1 = [0; abs(delta_hat)];
for j = 1:p1
    if delta1(j) == 1
        vlb1(j) = 1e-2;
        vub1(j) = 1;
    else 
        vlb1(j) = -inf;
        vub1(j) =  inf;
    end
end
vlb = [vlb1; vlb1];        
vub = [vub1; vub1];  
%% Slovling using interior point method
[X_sol, cost, exitflag, output, mu, grad, hessian] = fmincon(@(x)(costFunctionCox(x(1:p1),x(p1+1:2*p1))), [theta_0;u_0], A, b, Aeq, beq, vlb, vub, [], options);
theta = X_sol(1:p1);
u = X_sol(p1+1: 2*p1);

theta_hat = X_sol(2:p1); % coefficient
theta_hat0 = X_sol(1);   % intercept
u_hat0 = X_sol(p1+1);

%% Outpunt Cost and theta
% txt1 = importdata(['D:\E\博士\R_程序\SVM\Data\RTCGA\TCGA_pro_outcome_TN_log_comp_UNtest.txt']);  
fprintf('Cost at theta found by fmincon: %e', cost);
fprintf('\n')
para = [lambda_opt, cost]
%% Predict
% format long;
% txt1 = importdata(['Data_test/1.txt']);  
% test_data = txt1.data;
% [k,l] = size(test_data);
% X_test = test_data(1:k,2:l);   
% % y_test = test_data(1:k,1); 
% y_test = 2*(test_data(1:k,1)+1)-3; 
% [y_true, y_pre, AUC, RefereceResult] = PredictCox(X_test, y_test, theta); 
% disp(RefereceResult) 
% P_test = [y_pre, y_true];
% para = [lambda_opt, AUC, cost]
%% Save
csvwrite([path,'Result/Result',i,'/theta_hat',i,'.csv'],theta_hat);
csvwrite([path,'Result/Result',i,'/theta',i,'.csv'],theta);
csvwrite([path,'Result/Result',i,'/para',i,'.csv'],para);

