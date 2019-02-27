%
%AUTHOR: JUAN DIAZ AMORIN

function [theta] = deanoxsbm(Y)
[n,s] = size(Y);
Z = zeros(n,n+s);

% Objective function of the BCC model: min f'x = (0*lambda - sum(s/y));
% x = [lambda(nx1) ; s(mx1)]' 
% f is specific by each n DMU

% Constraints
Aeq = [Y' -eye(s)];   % (M+1) equality constraints
lb = zeros(n+m,1);                          % N + M inequalities (lower bounds for lambdas, s)

for j=1:n
    f = [zeros(1,n) (1+(1/s).*Y(j,:).^(-1))];
    beq = [Y(j,:)'];
    z = linprog(f,[],[],Aeq,beq,lb);
    Z(j,:) = z';
end
theta = Z(:,n+1:end);

end
