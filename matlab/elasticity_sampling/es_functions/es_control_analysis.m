function control = es_control_analysis(RSp, RJp, RSpp, RJpp, c, v, u, ind_ext)

% ES_CONTROL_ANALYSIS - Compute various control and response coefficients
%
% control = es_control_analysis(RSp, RJp, RSpp, RJpp, c, v, u, ind_ext)
%
% Given parameter response matrices RSp and RJp and tensors RSpp and RJpp, 
% compute the corresponding matrices and tensors for enzymes and external 
% metabolites separately and compute all scaled response coefficients
% 
% Distinguish between response coefficients for enzymes (u) and external 
% metabolites s = c(ind_ext), while p denotes all parameters p = [u; s]

nr = size(RJp,1);

% ---------------------------------------------------
% set scaled response coefficients

control.RSp_un  = RSp  ;
control.RJp_un  = RJp  ;
control.RSpp_un = RSpp ;
control.RJpp_un = RJpp ;

control.RSu_un  = RSp(:,1:nr);
control.RSs_un  = RSp(:,nr+1:end);
control.RSuu_un = RSpp(:,1:nr,1:nr);
control.RSus_un = RSpp(:,1:nr,nr+1:end);
control.RSss_un = RSpp(:,nr+1:end,nr+1:end);

control.RJu_un  = RJp(:,1:nr);
control.RJs_un  = RJp(:,nr+1:end);
control.RJuu_un = RJpp(:,1:nr,1:nr);
control.RJus_un = RJpp(:,1:nr,nr+1:end);
control.RJss_un = RJpp(:,nr+1:end,nr+1:end);


% ---------------------------------------------------
% compute scaled response coefficients

p = [u; c(ind_ext)];

[RSp_sc, RJp_sc, RSpp_sc, RJpp_sc] = norm_response_coefficients(c, v, p, RSp, RJp, RSpp, RJpp);

control.RSp_sc  = RSp_sc  ;
control.RJp_sc  = RJp_sc  ;
control.RSpp_sc = RSpp_sc ;
control.RJpp_sc = RJpp_sc ;

control.RSu_sc  = RSp_sc(:,1:nr);
control.RSs_sc  = RSp_sc(:,nr+1:end);
control.RSuu_sc = RSpp_sc(:,1:nr,1:nr);
control.RSus_sc = RSpp_sc(:,1:nr,nr+1:end);
control.RSss_sc = RSpp_sc(:,nr+1:end,nr+1:end);

control.RJu_sc  = RJp_sc(:,1:nr);
control.RJs_sc  = RJp_sc(:,nr+1:end);
control.RJuu_sc = RJpp_sc(:,1:nr,1:nr);
control.RJus_sc = RJpp_sc(:,1:nr,nr+1:end);
control.RJss_sc = RJpp_sc(:,nr+1:end,nr+1:end);
