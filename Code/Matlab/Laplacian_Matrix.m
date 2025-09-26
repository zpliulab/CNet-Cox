% Non.NormalizedLaplacianMatrix = function(adj){
%   diag(adj) <- 0
%   deg <- apply(adj,1,sum)
%   D = diag(deg)
%   L = D - adj             % The most common L matrix
%   return(L)
% }

function [L, L1] = Laplacian_Matrix(p,adj) % Non Normalized

for i = 1:p
    adj(i,i) = 0;
end

deg = sum(adj,2);           % Sum of elements in each row of the matrix
D = sparse(p,p);
for i = 1:p
    D(i,i) = deg(i);
%     D(i,i) = sqrt(deg(i));
end
L = D - adj;
L = sparse(L);

% D = diag(sum(adj,2)); 
D1 = sqrt(D); 
[dimd,~] = size(D1);
for i=1:dimd
    if ( D1(i,i)~= 0 )
        D1(i,i)=1/D1(i,i);
    end
end
L1 = D1*L*D1;
L1 = sparse(L1);
return
