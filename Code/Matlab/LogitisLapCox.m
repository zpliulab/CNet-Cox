function [theta_hat, theta_0] = LogitisLapCox (X,ytime,ystatus,R_matrix,L,lambda1,lambda2)
% global p 
%% To get the theta which minimize the cost function -sum_i{1-y_i(theta'*x_i)}
%% +lambda1 |theta|_1+lambda2 theta' A theta
%% The input y is a n by 1 matrix
%% Input X is G X n gene expression matrix
%% L is graphical Laplacian matrix, lambda1 and lambda2 are parameters in the cost function that are known
%% The output theta is the optimal solution of the covex optimization problem
%% The dimension of theta: G by 1
[dimn,~] = size(ytime); % get the number of n and set it to dimn
[dimG,m] = size(X); % get the number of G and set it to dimG

%%  To use cvx to solve the convex problem
cvx_begin
    variables theta(dimG+1);
    X_trainhat = [ones(dimn,1), X'];   % 605   545  
    z = X_trainhat * theta;     % 605   1
% It must vector*matrix. a*A=b. a and b is 1*3, A is 3*3
% Here z is 3*1, so it needs 2 times T. one for vector, and one for sum.
%     - ystatus_train * (X' * coef - log( sum( exp((X'*coef)' .* R_matrix_train), 2)' )')
%     minimize(  sum( log_sum_exp( z'* R_matrix') ) )  % right
% Here we dont use sum() function, we use ystatus vextor multiply vector
% Here can not use log_sum_exp, it is not rom 1 to N, add all of them,
% Here can use exp(z'*R) for sum in risk set and then log the sum
    minimize(  -ystatus' * ( z - log_sum_exp( z'* R_matrix') )+ lambda1*norm(theta(2:dimG+1),1)+lambda2*theta(2:dimG+1)'*L*theta(2:dimG+1) ); %% wrong
%     minimize(  -ystatus' * ( z - log( exp(z') * R_matrix')' )+ lambda1*norm(theta(2:dimG+1),1)+lambda2*theta(2:dimG+1)'*L*theta(2:dimG+1) );    
%     minimize(  -ystatus' * ( z - log_sum_exp( z'* R_matrix') ) ...
%     + lambda1*norm(theta(2:dimG+1),1)+lambda2*theta(2:dimG+1)'*L*theta(2:dimG+1) );
%%
cvx_end

theta_hat = theta(2:dimG+1);
theta_0 = theta(1);

return
