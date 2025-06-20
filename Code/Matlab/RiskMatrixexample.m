%% 2023.4.14 LLY created
% This file is to validate the effective of the objection
% According to an example list in my note
% There are three people dead at time 1, time 3 and time 5.
% The s is the loss function value by manul, and the last line is by coding
% The key point is the R_matirx and the sum(vector,2). the vector multiply
% each row of the matrix

%% clear
clc; clear; close all; format long;


%% 
% theta = np.array([1,2,3]).reshape(1,3)
% R = np.array([[1,1,1],[0,1,1],[0,0,1]])
% print(np.sum(theta*R,axis=1))

%% vector and matrix
theta = [1 2 3];
R = [1 1 1; 0 1 1; 0 0 1];
size(theta)    % 1*3
size(R)     % 3*3
% theta*R
% theta.*R
% sum(theta*R, 1)
% python 中的*是对应位置相乘（矩阵用np.dot（）），
% matlab中的*是矩阵乘法，对应位置相乘是 .*
theta .* R
sum(theta .* R, 2)
theta * R'
sum(R .* theta, 2)

sum(exp(theta) .* R, 2)

sum(theta * R, 2)

delta = [1 1 1];
size(delta);
delta * (theta - sum(theta .* R, 2))


c = zeros(3);
for i=1:size(R,1)
  c(i,:) = R(i,:) .* theta;
end
%% data matrix and coefficients
X = [1 2 3; 4 5 6; 7 8 9]';
X1 = X(:,1);
X2 = X(:,2);
X3 = X(:,3);
size(X1)    % 3*1
coef = [0.3 0.6 0.9]'
size(coef)
theta1 = X1'*coef;
theta2 = X2'*coef;
theta3 = X3'*coef;

%% matrix multiply vector
theta = [theta1; theta2; theta3]
ystatus_train = [1 1 1];
size(X')
size(coef)
X' * coef
%% compute the example manual
s1 = ystatus_train(1,1) * (theta1 - log(exp(theta1) + exp(theta2) + exp(theta3)));
s2 = ystatus_train(1,2) * (theta2 - log(exp(theta2) + exp(theta3)));
s3 = ystatus_train(1,3) * (theta3 - log(exp(theta3)));
[s1 s2 s3]
[log(exp(theta1) + exp(theta2) + exp(theta3)), log(exp(theta2) + exp(theta3)), log(exp(theta3))]
[exp(theta1) + exp(theta2) + exp(theta3), exp(theta2) + exp(theta3), exp(theta3)]
s = s1 + s2 + s3
%% survival time and risk set
ystatus_train = [1 1 1];
ytime_train = [1 2 3];
N = size(ytime_train);
R_matrix_train = sparse(N(1,2), N(1,2));
for i = 1:N(1,2)
    for j = 1:N(1,2)
%     for j = i:N(1,2)
        if ytime_train(1,j) >= ytime_train(1,i)
            R_matrix_train(i,j) = 1;
        end
    end
end
%% visualization
spy(R_matrix_train)
%% compute the example coding, we need to maximun the log partical likelihood function
ystatus_train * (X' * coef - log( sum( exp((X'*coef)' .* R_matrix_train), 2)' )')

%% The objective function need to be negtivate, we need to minimize
- ystatus_train * (X' * coef - log( sum( exp((X'*coef)' .* R_matrix_train), 2)' )')


%% Test matrix * matrix
aa = [1 2 3; 4 5 6];
bb = R;
aa.*bb

