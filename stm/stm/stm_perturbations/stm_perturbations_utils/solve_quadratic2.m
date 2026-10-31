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
% 0   = x' A   + l' B + l' (x' D )
% 0   = ( A + l' D ) x + B' l
% x = - inv( A + l o D) B' l
%
% inserting it into (2) yields
% 0 = - B ( inv( A + l o D) B' l) + 1/2 D * (( inv( A + l o D) B' l) o ( inv( A + l o D) B' l)) - c
% which has to be solved for l
%
% for the starting guess (omit quadratic term with D):
%
% min = 1/2 x' A x + l' B x
% => x = - inv(A) B' l
% insert in (2) in the form "B*x = c"
% yields - B inv(A) B' l  = c    so guess  l = -inv(B inv(A) B') c
%
% Example:
%   A = eye(3);
%   B = [1 1 1; 0 0 1];
%   c = [1 2]';
%   D(1,:,:) = eye(3);
%   D(2,:,:) = eye(3);
%   x = solve_quadratic2(A,c,B,D)

eval(default('seed','0'));

% use a fixed seed in this function only (the state of the generator is restored below)
rng_state = rng;
rng(seed);

l_guess = - pinv(B*pinv(A)*B')*c;
l_guess = l_guess + randn(size(l_guess));

l       =   fsolve(@(l) solve_quadratic2_fct(l,A,B,c,D), l_guess, optimoptions('fsolve','Display','off'));
x       = - inv( A + squeeze(tensor_product(l',D,2,1)) ) * B' * l;

rng(rng_state);

cpred   = B * x + 1/2 * tensor_product(tensor_product(D,x,2,1),x,2,1);

mismatch = abs( (cpred - c) ./ (cpred + c) );

if max(mismatch)>10^-10
  warning('solve_quadratic2: constraint mismatch %g', max(mismatch));
end
