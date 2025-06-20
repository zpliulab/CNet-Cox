function [lambda_opt, alpha_opt, r] = cvCox(X,ytime,y, R_matrix, L,alpha,lambda,nfold)
%% The function for logistic regression network regularization applications  
%% The function is to use cross validation (CV) to get the optimal value of
%% parameter lambda and alpha through the training and testing procedures.
%% nfold is the number of folds in CV, which could be selected randomly
%% For Lasso case, alpha is 1, which is included in the code instruction.
% X = X_train;
% ytime = ytime_train;
% y = ystatus_train;
% L = L1;

X = X';
[dimG,~]=size(X); % get the dimension G 
[dimn,~]=size(y); % get the dimension n
%% Get the number of candidate lambda's
[m,n]=size(lambda);
if m>n
    sizelam=m;
else
    sizelam=n;
end
%% Get the length of alpha. for lasso: alpha=1
[m,n]=size(alpha);
if m>n
    sizealpha=m;
else
    sizealpha=n;
end
%% To divide X and y into several folds respectively
for i=1:(nfold-1)
  fnw=['x',int2str(i),'=X(:,(floor(dimn/nfold)*',int2str(i),'-floor(dimn/nfold)+1):(floor(dimn/nfold)*',int2str(i),'));'];
  eval(fnw);
  fnw=['y',int2str(i),'=y((floor(dimn/nfold)*',int2str(i),'-floor(dimn/nfold)+1):(floor(dimn/nfold)*',int2str(i),'),1);'];
  eval(fnw);
%   % LLY added 2023.4.14
  fnw=['ytime',int2str(i),'=ytime((floor(dimn/nfold)*',int2str(i),'-floor(dimn/nfold)+1):(floor(dimn/nfold)*',int2str(i),'),1);'];
  eval(fnw);
% calculate the Risk matrix under the selected fold ytime
  fnw=['R_matrix',int2str(i),'= RiskMatrix(ytime',int2str(i),');'];
  eval(fnw);
i
end
 
%% The last one is the columns that left 
fnw=['x',int2str(nfold),'=X(:,(floor(dimn/nfold)*',int2str(nfold-1),'+1):dimn);'];
eval(fnw);
fnw=['y',int2str(nfold),'=y((floor(dimn/nfold)*',int2str(nfold-1),'+1):dimn);'];
eval(fnw);
% 2023.4.14 LLy added
fnw=['ytime',int2str(nfold),'=ytime((floor(dimn/nfold)*',int2str(nfold-1),'+1):dimn);'];
eval(fnw);
fnw=['R_matrix',int2str(nfold),'= RiskMatrix(ytime',int2str(nfold),');'];
eval(fnw);
%% To set the head and rear xnot first.
xnot1=X(:,floor(dimn/nfold)+1:dimn);
fnw=['xnot',int2str(nfold),'=[X(:,1:(floor(dimn/nfold)*',int2str(nfold-1),'))];'];
eval(fnw);

%% Set the rest xnot's
for i=2:(nfold-1)
fnw=['xnot',int2str(i),'=[X(:,1:(floor(dimn/nfold))*(',int2str(i),'-1))  X(:, (floor(dimn/nfold)*',int2str(i),'+1):dimn)];'];
eval(fnw);
end

%% Set ynot in a similar way -- old
ynot1=y(floor(dimn/nfold)+1:dimn,1);
fnw=['ynot',int2str(nfold),'=[y(1:(floor(dimn/nfold)*',int2str(nfold-1),'),1)];'];
eval(fnw);
for i=2:(nfold-1)
fnw=['ynot',int2str(i),'=[y(1:(floor(dimn/nfold))*(',int2str(i),'-1)); y((floor(dimn/nfold)*i+1):dimn, 1)];'];
eval(fnw);
end

%% calcus the y_time
% 2023.4.14 LLY add, just for y time splitting
ytimenot1 = ytime(floor(dimn/nfold)+1:dimn,1);
R_matrixnot1 = RiskMatrix(ytimenot1);
fnw=['ytimenot',int2str(nfold),'=[ytime(1:(floor(dimn/nfold)*',int2str(nfold-1),'),1)];'];
eval(fnw);
fnw=['R_matrixnot',int2str(nfold),'= RiskMatrix(ytimenot',int2str(nfold),');'];
eval(fnw);
for i=2:(nfold-1)
fnw=['ytimenot',int2str(i),'=[ytime(1:(floor(dimn/nfold))*(',int2str(i),'-1)); y((floor(dimn/nfold)*i+1):dimn, 1)];'];
eval(fnw);
fnw=['R_matrixnot',int2str(i),'= RiskMatrix(ytimenot',int2str(i),');'];
eval(fnw);
i
end
%% The following is the cross validation to choose optimal alpha and lambda
for k=1:sizealpha
alphacan = alpha(k);    
% alphacan = alpha(1);  
for i=1:nfold
%     i = 1;
 for j=1:sizelam
%      j = 1;
% 2023.4.14 LLY added 
% eval函数作用简单来说就是可以把字符串当作命令来执行。即将字符串自动识别并转化为matlab命令
     fnw=['[Theta,Theta0]=LogitisLapCox(xnot',int2str(i),',ytimenot',int2str(i),',ynot',int2str(i),', R_matrixnot',int2str(i),', L,lambda(',int2str(j),')*alphacan,lambda(',int2str(j),')*(1-alphacan));'];
     eval(fnw);
     %% r is the matrix that contains the residual=|| yi-(1/(1+exp(-theta'*Xi))>0.5) ||_2
%     fnw=['r(',int2str(i),',',int2str(j),')=norm((transpose(y',int2str(i),')-(1./(1+exp(-transpose(Theta)*x',int2str(i),'))>.5)),2);'];
%     fnw=['r(',int2str(i),',',int2str(j),')=norm((transpose(y',int2str(i),')-sign(transpose(Theta)*x',int2str(i),'+ Theta0)),2);'];
%  y1'*( transpose(x1)*Theta+Theta0 -transpose(log( exp(transpose(Theta)*x1+Theta0) * transpose(R_matrix1) )))
%  transpose(ystatus_train)*((X_train*Theta + Theta0) - transpose(log( exp(transpose(Theta)*transpose(X_train)+Theta0) *transpose(R_matrix) )))
    fnw=['r(',int2str(i),',',int2str(j),')=transpose(y)*((transpose(X)*Theta + Theta0) - transpose(log( exp(transpose(Theta)*X+Theta0) *transpose(R_matrix) ))) - transpose(y',int2str(i),') *(transpose(x',int2str(i),')*Theta+Theta0 - transpose( log( exp(transpose(Theta)*x',int2str(i),'+ Theta0) * transpose(R_matrix',int2str(i),')) ));'];
    eval(fnw); 
 end
end

%% to get the minimal residual and the corresponding index, which is the index for optimal lambda
rmin(k)= max(mean(r,1)); % first set the maximal mean to rmax
optindex=1;  % in case, the optimal lambda is the first one
for j=1:sizelam
 if rmin>mean(r(:,j))
     rmin(k)=mean(r(:,j)); 
     optindex=j;  % if mean is less, take the index and reset the minimal residual
 end
end
%% Get the optimal lambda
lambdaopt(k)=lambda(optindex);  % For every alpha, there will be an optimal lambda 
end
[err Ind]=min(rmin); 
alpha_opt=alpha(Ind);
lambda_opt=lambdaopt(Ind);

    