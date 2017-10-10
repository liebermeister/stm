function x = solve_quadratic(A,c,B)

% x = solve_quadratic(A,B,c)
% solve min = x'*A*x with constraint B*x = c
% with symmetric matrix A

K      = null(B);
pinv_B = pinv(full(B));
Id     = eye(size(A,1));
x      = (Id - K*inv(K'*A*K)*K'*A) * pinv_B * c;

% proof: set x = pinv(B)*c + K * w
% with B*K = 0 to satisfy constraint
% insert into quadratic form and solve for the minimum
% test
% A = eye(3);
% B = [1 1 1; 0 0 1];
% c = [1 2]';
% x = solve_quadratic(A,c,B);
% B * x;
