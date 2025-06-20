function J = costFunction13(theta, u)
global  L1 X_train  ystatus_train  R_matrix  lambda_opt alpha_opt m p1

%%
% m = length(y); % number of training examples
ystatus = ystatus_train;
X_trainhat = [ones(m,1), X_train];
z = X_trainhat * theta;
% hx = 1 ./ (1 + exp(-z));

u_hat = u(2:p1);
R1 = sum(u_hat);
% R1 = sum(abs(theta));
theta_hat = theta(2:p1);
R2 = theta_hat'*L1*theta_hat/2;
%%
% J = sum([- y_train' * log(hx) - (1 -  y_train)' * log(1 - hx)]) + lambda_opt*alpha_opt*R1 + lambda_opt*(1-alpha_opt)*R2;
% J = sum( max(0, 1-y_train'*z) ) + lambda_opt*alpha_opt*R1 + lambda_opt*(1-alpha_opt)*R2;
J = - ystatus' * ( z - log_sum_exp( z'* R_matrix') ) + lambda_opt*alpha_opt*R1 + lambda_opt*(1-alpha_opt)*R2;

return