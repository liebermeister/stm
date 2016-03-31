function r = es_simulate_perturbation(N,W,ind_ext,u_ratio,s_ratio,es_result,expansion,expansion_order)

%ES_SIMULATE_PERTURBATION - Predict the effects of a given differential profile 
%
% r = es_simulate_perturbation(N,W,ind_ext,u_ratio,s_ratio,es_result,expansion,expansion_order)
%
% Compute the predicted effects of a given differential profile 
% of enzyme expression (u) and external metabolites (s)
% on balanced concentrations c and fluxes j
% based on an expansion using the response coefficients in "es_result"
% 
% N, W, ind_ext: structural information about the network
% s_ratio:       external concentration ratio vector
% u_ratio:       expression ratio vector
% es_result:     output from a previous elasticity sampling (function 'sample_model')
% expansion:     {'logarithmic','non-logarithmic'}
%                -> depending on the expansion type, the output refers to logarithmic 
%                   or non-logarithmic values; 
% expansion_order          1 or 2 expansion_order of the expansion
%
% output r contains the ratios etc.

p_ratio = [u_ratio; s_ratio];

[nm,nr] = size(N);
nsample = size(u_ratio,2);
c       = es_result.c;
v       = es_result.v;

switch expansion,
  
  case 'logarithmic',

% fits better the multiplicative nature of kinetic models
% but doesn't preserve stationarity condition of fluxes
    
    RSpsc   = es_result.control.RSp_sc;
    RJpsc   = es_result.control.RJp_sc;
    
    switch expansion_order,
      
      case 1
        
        d_c = squeeze(RSpsc(:,:)) * log(p_ratio);
        d_v = squeeze(RJpsc(:,:)) * log(p_ratio);
        
      case 2
        
        RSppsc  = es_result.control.RSpp_sc;
        RJppsc  = es_result.control.RJpp_sc;
        
        for i2 = 1:nsample          
          dlogS2(:,i2) =  tensor_product(tensor_product(...
              squeeze(RSppsc(:,:,:)),...
              log(p_ratio(:,i2)),2,1),log(p_ratio(:,i2)),2,1);
          dlogJ2(:,i2) =  tensor_product(tensor_product(...          
              squeeze(RJppsc(:,:,:)),...
              log(p_ratio(:,i2)),2,1),log(p_ratio(:,i2)),2,1);
        end
        
        d_c = squeeze(RSpsc(:,:)) * log(p_ratio) + 0.5 * dlogS2;
        d_v = squeeze(RJpsc(:,:)) * log(p_ratio) + 0.5 * dlogJ2;
    end

    d_c(ind_ext) = log(s_ratio);
    
  case 'non-logarithmic',
    
% preserves stationarity condition of fluxes
% but doen't fit the multiplicative nature of kinetic models
    
    u = es_result.u;
    s = es_result.c(ind_ext);
    p = [u; s];
    delta_p = p_ratio .* p - p;
    
    RSp   = es_result.control.RSp;
    RJp   = es_result.control.RJp;
    
    switch expansion_order
      
      case 1
        
        d_c = squeeze(RSp(:,:)) * log(p_ratio);
        d_v = squeeze(RJp(:,:)) * log(p_ratio);
        
      case 2
        
        RSpp  = es_result.control.RSpp;
        RJpp  = es_result.control.RJpp;
        
        for i2 = 1:nsample
          delta_S2(:,i2) =  tensor_product(tensor_product(...
              squeeze(RSpp(:,:,:)),...
              delta_p(:,i2),2,1),delta_p(:,i2),2,1);
          delta_J2(:,i2) =  tensor_product(tensor_product(...
              squeeze(RJpp(:,:,:)),...
              delta_p(:,i2),2,1),delta_p(:,i2),2,1);
        end
        
        d_c = squeeze(RSp(:,:)) * log(p_ratio) + 0.5 * delta_S2;
        d_v = squeeze(RJp(:,:)) * log(p_ratio) + 0.5 * delta_J2;
        
    end
    
    d_c(ind_ext) = s_ratio .* s - s;
    
end


% ----------------------------------------------------
% collect the results in the output data structure r

r.c_init = es_result.c;
r.v_init = es_result.v;

switch expansion, 
  
  case ('non-logarithmic'), 

    % exclude negative concentrations
    d_c(r.c_init + d_c<0) = -r.c_init(r.c_init + d_c<0);
    r.c_ratio  = 1 + d_c./r.c_init;
    r.c_new    = r.c_init + d_c;
    r.v_ratio  = 1 + d_v./r.v_init;
    r.v_new    = r.v_init + d_v;
  
  case ('logarithmic'),

    r.c_ratio  = exp(d_c);
    r.c_new    = r.c_ratio .* r.c_init;
    % ensure stationarity: 
    K = null(N(find(ind_ext==0),:)); % kernel matrix
    dum = exp(d_v); dum(isnan(dum))=0; 
    d_v_nonlog = K * pinv(K) * [r.v_init .* dum];
    r.v_ratio  = 1 + d_v_nonlog./r.v_init;
    r.v_new    = r.v_init + d_v_nonlog;

end

r.delta_c  = r.c_new - r.c_init;
r.delta_v  = r.v_new - r.v_init;
