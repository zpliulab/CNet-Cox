function [theta_hat, theta_0] = LogitisLapCox (X,ytime,ystatus,R_matrix,L,lambda1,lambda2)
% global p 
%% To get the theta which minimize the cost function -sum_i{1-y_i(theta'*x_i)}
%% +lambda1 |theta|_1+lambda2 theta' A theta
%% The input y is a n by 1 matrix
%% Input X is G X n gene expression matrix
%% L is graphical Laplacian matrix, lambda1 and lambda2 are parameters in the cost function that are known
%% The output theta is the optimal solution of the covex optimization problem
%% The dimension of theta: G by 1
X = xnot1;
ytime = ytimenot1;
ystatus = ynot1;
R_matrix = R_matrix1;

[dimn,~] = size(ytime); % get the number of n and set it to dimn
[dimG,m] = size(X); % get the number of G and set it to dimG

%%  To use cvx to solve the convex problem

X_trainhat = [ones(dimn,1), X'];   % 605   545  
    
cvx_begin
    variables theta(dimG+1);
%     variables theta(dimn);

%% 测试 max 函数，可行
    z = X_trainhat * theta;     % 605   1
%% 加上惩罚项
% It must vector*matrix. a*A=b. a and b is 1*3, A is 3*3
% Here z is 3*1, so it needs 2 times T. one for vector, and one for sum.
%%
% bb = sum(z' .* R_matrix, 2);   %% Matrix dimensions must agree.
%%
%     variables z(dimn)
%     expression b(dimn,1)
%     b = zeros(dimn, 1);
%     for i = 1:+
%         dimn
%         for j = 1:dimn
%             b(i) = b(i) + R_matrix(i,j) * z(j); % 无法从 cvx 转换为 double
%         end
%     end
%%
%     ystatus_train * (X' * coef - log( sum( exp((X'*coef)') .* R_matrix_train, 2)' )')
%     minimize(  sum( exp(z') )  )    % right
%     minimize(  sum( exp(z') * R_matrix))     % right 向量求和
%     minimize(  sum( (exp(z))'* R_matrix))     % right
%     minimize(  sum( log_sum_exp( X_trainhat * theta) ) )  % right
%     minimize(  sum( log_sum_exp( b') ) )  % wrong
    minimize(  sum( log_sum_exp( z'* R_matrix') ) )  % wrong
%     minimize( (exp(z))' .* R_matrix )    % wrong
%     minimize( sum(theta' .* R_matrix, 2) )    % wrong
%     sum( log_sum_exp( ((X_trainhat * theta) .* R_matrix), 2) )    % wrong
%     minimize(  sum(sum( (exp(z))' .* R_matrix, 2) ))    % wrong
%     minimize(  sum( sum( (exp(z))' .* R_matrix, 2) )  )    % right
%   minimize( sum( max( zeros(m,1), ones(m,1)-ystatus.*z ) )) 
%     minimize(  sum( sum( (exp(z_hat))' .* R_matrix, 2)' )  )    % right
%     + lambda1*norm(theta(2:dimG+1),1)+lambda2*theta(2:dimG+1)'*L*theta(2:dimG+1) );
%%
cvx_end

theta_hat = theta(2:dimG+1);
theta_0 = theta(1);

return
%% the same with the following cycle. Just one element

%     theta0([1:(dimG+1)],1) = -1; 
%     z_hat = X_trainhat * theta0;     % 605   1 

% aa = sum(z_hat' .* R_matrix(1,:), 2); 
% a = zeros(dimn, 1);
% % for i = 1:dimn
%     i = 1;
%     for j = 1:dimn
% %         j = 1
%         a(i) = a(i) + z_hat(j) * R_matrix(i,j);
%     end
%% the same with the following cycle. The whole mattrix

%     theta0([1:(dimG+1)],1) = -1; 
%     z_hat = X_trainhat * theta0;     % 605   1 

% bb = sum(z_hat' .* R_matrix, 2); 
% b = zeros(dimn, 1);
% for i = 1:dimn
% %     i = 1;
%     for j = 1:dimn
% %         j = 1
%         b(i) = b(i) + z_hat(j) * R_matrix(i,j);
%     end
% end

%% .* matrix  is rqualent to * matrix's transport
% cc = z_hat' * R_matrix';     % 1*605
% cc  
% sum(cc)
% sum(cc,2)    % 3.8527e+05
% dd = sum(z_hat' .* R_matrix, 2)'    % 1*605
% sum(dd,2)    % 3.8527e+05
