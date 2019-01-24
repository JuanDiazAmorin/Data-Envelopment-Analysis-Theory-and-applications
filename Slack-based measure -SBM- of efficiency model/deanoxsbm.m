%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEANOX SBM
% 
% Implements radial DEA without inputs, slacks-based measure efficiency
% See pdf. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data envelopment analysis.
% Only outputs, no inputs. 
% Y: nxs vector of outputs
% theta: nxm matrix of slacks (one line for each DMU)

function [theta] = deanoxsbm(Y)

[n,m] = size(Y);

Z = zeros(n,n+m);

% Objective function of the BCC model: min f'x = (0*lambda - sum(s/y));
% x = [lambda(nx1) ; s(mx1)]' 
% f is specific by each n DMU

% Constraints
Aeq = [Y' -eye(m); ones(1,n) zeros(1,m)];   % (M+1) equality constraints
lb = zeros(n+m,1);                          % N + M inequalities (lower bounds for lambdas, s)

for j=1:n
    f = [zeros(1,n) -(1/m).*Y(j,:).^(-1)];
    beq = [Y(j,:)'; 1];
    z = linprog(f,[],[],Aeq,beq,lb);
    Z(j,:) = z';
end
theta = Z(:,n+1:end);

end







