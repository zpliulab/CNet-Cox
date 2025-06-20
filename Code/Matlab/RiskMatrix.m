%% 2023.4.14 LLY created
% This file is to validate the effective of the objection
% According to an example list in my note
% There are three people dead at time 1, time 3 and time 5.
% The s is the loss function value by manul, and the last line is by coding
% The key point is the R_matirx and the sum(vector,2). the vector multiply
% each row of the matrix

function [R_matrix] = RiskMatrix(ytime)
% survival time and risk set
% ytime_train = ytime_train(1:4);
ytime_train = ytime;
N = size(ytime_train);
[N,~] = size(ytime_train);
R_matrix_train = sparse(N, N);
for i = 1:N
    for j = 1:N
        if ytime_train(j) >= ytime_train(i)
            R_matrix_train(i,j) = 1;
        end
    end
end
R_matrix = R_matrix_train;
% spy(R_matrix_train);
return

