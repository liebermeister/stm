function [u_ratio_pred, s_ratio_pred, pred_effect] = es_discover_perturbation(N,W,ind_ext,c_ratio,v_ratio,es_result,expansion,expansion_order,u_sigma_prior,s_known_sigma_prior,s_unknown_sigma_prior)

% [u_ratio_pred, s_ratio_pred, pred_effect] = es_discover_perturbation(N,W,ind_ext,c_ratio,v_ratio,es_result,expansion,expansion_order,u_sigma_prior,s_known_sigma_prior,s_unknown_sigma_prior)
%
% Compute the most likely perturbation of enzyme expression (u) and external metabolites (s)
% that would give rise to a given differential profile of balanced concentrations c and fluxes j
% 
% This is, more or less, an inversion of 'es_simulate_perturbation'
%
% N, W, ind_ext: structural information about the network
% c_ratio,v_ratio vectors containing perturbations (no matrices for several samples allowed!)
% es_result:     output from a previous elasticity sampling (function 'sample_model')
% expansion:     {'logarithmic','non-logarithmic'}
%                -> depending on the expansion type, the output refers to logarithmic 
%                   or non-logarithmic values; 
% expansion_order:       1 or 2 for first- or second-expansion_order expansion
%
% output:
%
% s_ratio:       external concentration ratio vector
% u_ratio:       expression ratio vector

sigma_decrease_factor = 0.8;

[nm,nr] = size(N);
nsample = size(c_ratio,2);
c       = es_result.c;
v       = es_result.v;
np      = nr+length(ind_ext);

if exist('u_sigma_prior','var'),
  dum = isfinite(c_ratio);
  dum = dum(ind_ext);
  s_sigma_prior = s_known_sigma_prior * dum + s_unknown_sigma_prior * (1-dum);
  C_prior_inv = diag( 1./[u_sigma_prior.^2 * ones(nr,1); s_sigma_prior.^2]);
else
  C_prior_inv = eye(np);
end

switch expansion 
  
  case 'logarithmic',

    % fits better the multiplicative nature of kinetic models
    % but doesn't preserve stationarity condition of fluxes
    
    RSpsc   = es_result.control.RSp_sc;
    RJpsc   = es_result.control.RJp_sc;

    % merge the known fold changes to a vector y and compute the response coeffs  
    
    d_log_y = [log(c_ratio); log(v_ratio)];
    RYpsc   = [RSpsc; RJpsc]; 
    ind_known = find(isfinite(d_log_y));
    d_log_y = d_log_y(ind_known);
    
    switch expansion_order,
      
      case 1,
        d_log_p = inv(RYpsc' * RYpsc + C_prior_inv) * RYpsc' * d_log_y;
        
      case 2,
        RSppsc  = es_result.control.RSpp_sc;
        RJppsc  = es_result.control.RJpp_sc;
        RYppsc(1:nm,:,:)      = RSppsc;
        RYppsc(nm+(1:nr),:,:) = RJppsc;
        RYpsc   = RYpsc(ind_known,:);
        RYppsc  = RYppsc(ind_known,:,:);        
        d_log_p = solve_quadratic2(C_prior_inv,d_log_y,RYpsc,RYppsc);
    
    end
    
    p_ratio = exp(d_log_p);
    
  case 'non-logarithmic',
    
    % preserves stationarity condition of fluxes
    % but doen't fit the multiplicative nature of kinetic models

    delta_c = (c_ratio - 1) .* c;
    delta_v = (v_ratio - 1) .* v;
    u = es_result.u;
    s = es_result.c(ind_ext);
    p = [u; s];
    
    RSp   = es_result.control.RSp;
    RJp   = es_result.control.RJp;

    % merge the known differences to a vector y and compute the response coeffs  
    
    delta_y   = [delta_c; delta_v];
    RYp       = [RSp; RJp]; 
    ind_known = find(isfinite(delta_y));
    delta_y   = delta_y(ind_known);
    RYp       = RYp(ind_known,:);
        
    switch expansion_order,
      case 1

        delta_p = inv( RYp' * RYp + C_prior_inv ) * RYp' * delta_y;

      case 2
        RSpp  = es_result.control.RSpp;
        RJpp  = es_result.control.RJpp;
        RYpp(1:nm,:,:)      = RSpp;
        RYpp(nm+(1:nr),:,:) = RJpp;
        RYpp    = RYpp(ind_known,:,:);
        delta_p = solve_quadratic2(C_prior_inv, delta_y, RYp, RYpp);
    
    end

    p_ratio = 1 + delta_p ./ p;

end

u_ratio_pred = p_ratio(1:nr);
s_ratio_pred = p_ratio(nr+1:end);

if sum(p_ratio<0), 
  display(sprintf(' Expansion leads to negative values. Decreasing the prior sigmas by a factor %f',sigma_decrease_factor)); 
  [u_ratio_pred, s_ratio_pred] = es_discover_perturbation(...
      N,W,ind_ext,c_ratio,v_ratio,es_result,expansion,...
      expansion_order,sigma_decrease_factor*u_sigma_prior,sigma_decrease_factor*s_known_sigma_prior,sigma_decrease_factor*s_unknown_sigma_prior);
end


if nargout >2,
  
  pred_effect = es_simulate_perturbation(...
      N,W,ind_ext,u_ratio_pred,s_ratio_pred,es_result,expansion,expansion_order);
  
end
