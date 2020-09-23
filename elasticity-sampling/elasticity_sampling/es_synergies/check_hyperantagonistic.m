function [success, solution] = check_hyperantagonistic(fa,fb,faa,fab,fbb)

% [success, solution] = check_hyperantagonistic(fa,fb,faa,fab,fbb)
%
% fitness function  f(a,b) = f0 + fa a + fb b + 1/2 (faa a^2  + 2 fab a b + fbb b^2)
% requirements: f(-a,0) < f(0,0); f(0,-b) < f(0,0); f(-a,-b) > f(-a,0); f(-a,-b) > f(0,-b); 

% test 
% fa = 1; fb = 1; faa=0; fab= 1; fbb = 0;  [success,solution] = check_hyperantagonistic(fa,fb,faa,fab,fbb)

solution = [nan; nan];

success = 0;
stop    = 0;

if faa < 0, 
  if fa < 0, a_range  = [2 * abs(fa/faa), inf]';
  else       a_range  = [0, inf]';
  end
else,
  if fa < 0, a_range  = [nan nan]'; stop = 1;
  else       a_range  = [0, 2 * fa/faa]';
  end
end

if fbb < 0, 
  if fb < 0, b_range  = [2 * abs(fb/fbb), inf]';
  else       b_range  = [0, inf]';
  end
else,
  if fb < 0, b_range  = [nan nan]'; stop = 1;
  else       b_range  = [0, 2 * fb/fbb]';
  end
end

if stop ==0,
  M = [1 0; 0 1; -1 0; 0 -1; faa, 2*fab; 2* fab, fbb];
  x = [a_range(1); b_range(1);-a_range(2);-b_range(2);fa;fb];
  relevant = find(isfinite(x));
  M = -M(relevant,:);
  a = -x(relevant)+10^-10;

[solution,goal,status] = my_linear_programming(M,a,ones(size(M,2),1));

solution(isnan(solution)) = 1;
if strcmp(status,'optimal'),
  if sum((M*solution <a)==0) ==0,
    display('Hyperantagonistic solution found');
    success = 1;
  else
    solution = nan *solution;
  end
else,
  solution = nan *solution;
end

end
