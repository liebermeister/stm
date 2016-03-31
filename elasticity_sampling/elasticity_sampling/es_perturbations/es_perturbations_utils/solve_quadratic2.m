function x = solve_quadratic2(A,c,B,D,seed)

% x = solve_quadratic2(A,c,B,D,seed)
% solve min = 1/2 x'*A*x with constraint (1)
% where B*x + 1/2 D (x o x ) = c     (2)
% with symmetric matrix A and symmetric tensor D
%
% with a langrangian multiplier vector l 
% we get 
%
% min = 1/2 x' A x + l' ( B x + 1/2 D (x o x) )
% 0   = x' A   + l' B + 1/2 l' (x' D )
% 0   = ( A + 1/2 l' D ) x + B' l
% x = - inv( A + 1/2 l o D) B' l 
%
% inserting it into (2) yields
% 0 = - B ( inv( A + 1/2 l o D) B' l) + 1/2 D * (( inv( A + 1/2 l o D) B' l) o ( inv( A + 1/2 l o D) B' l)) - c
% which has to be solved for l
%
% for the starting guess (omit quadratic term with D):
%
% min = 1/2 x' A x + l' B x
% => x = - inv(A) B' l
% insert in (2) in the form "B*x = c"
% yields - B inv(A) B' l  = c    so guess  l = -inv(B inv(A) B') c

eval(default('seed','0'));
randn('state',seed);

l_guess = - pinv(B*pinv(A)*B')*c;
l_guess = l_guess + randn(size(l_guess));

l       =   fsolve(@solve_quadratic2_fct,l_guess,[],A,B,c,D);
x       = - inv( A + 1/2 * squeeze(tensor_product(l',D,2,1)) ) * B' * l;

cpred   = B * x + 1/2 * tensor_product(tensor_product(D,x,2,1),x,2,1);

mismatch = abs( (cpred - c) ./ (cpred + c) );

if max(mismatch)>10^-10,
  max_mismatch = max(mismatch)
end

% Test

return

A = eye(3);
B = [1 1 1; 0 0 1];
c = [1 2]';
D(1,:,:) = eye(3);
D(2,:,:) = eye(3);

x = solve_quadratic2(A,c,B,D)
