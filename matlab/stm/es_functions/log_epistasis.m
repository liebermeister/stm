function epistasis = simple_epistasis(f_ref,f_single,f_double)

% LOG_EPISTASIS Simple epistasis measure
%
% epistasis = simple_epistasis(f_ref,f_single,f_double)
%
% Input 
% f_ref    wildtype fitness (number)
% f_single single perturbation fitnesses (vector)
% f_double double perturbation fitnesses (matrix)
%
% Output
% epistasis (matrix)

epsilon = 10^-8;

w_single = f_single/f_ref;
w_double = f_double/f_ref;
w_single(abs(w_single)<10^-8)= 0;
w_double(abs(w_double)<10^-8)= 0;

w_ratio = w_double ./ [w_single * w_single'];
w_ratio(w_double == [w_single * w_single']) = 1;
epistasis = w_ratio-1;