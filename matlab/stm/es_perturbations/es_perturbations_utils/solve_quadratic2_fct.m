function f = solve_quadratic2_fct(l,A,B,c,D)

xx = inv( A + 1/2 * squeeze(tensor_product(l',D,2,1)) ) * B' * l;

f  = - B * xx + 1/2 * tensor_product(tensor_product(D,xx,3,1),xx,2,1) - c;
