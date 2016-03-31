%METABOLIC_CORRELATIONS_GRAPHICS Graphics for metabolic correlations 
% 
% metabolic_correlations_graphics(network,R_list,C_list,corr_list,fignum)


function metabolic_correlations_graphics(network,R_list,C_list,corr_list,fignum)

eval(default('fignum','1:3'));

[nr,nm,nx] = network_numbers(network);
n_int      = sum(network.external==0);
nit        = size(R_list,1);

internal = find(network.external ==0);


% ------------------------------------------------------
% distribution of absolute response coefficients

R_log_list = log(abs(R_list));
R_log_mean = squeeze(mean(R_log_list,1));
R_log_std  = reshape( sqrt(mean( (reshape(R_log_list,nit,n_int*nr) - repmat(reshape(R_log_mean,1,n_int*nr),nit,1)).^2 )),n_int,nr);

R_sign_mean = squeeze(mean(sign(R_list),1));

max_value = max(max(exp(R_log_mean+R_log_std)));

figure(fignum(1)); clf
im_circle( exp(R_log_mean+R_log_std) .*(R_sign_mean>0),max_value,network.metabolites(internal),network.actions,[.5 .5 1]); hold on;
im_circle( exp(R_log_mean) .*(R_sign_mean>0),max_value,network.metabolites(internal),network.actions,[0 0 1]);  hold on;
im_circle( exp(R_log_mean+R_log_std) .*(R_sign_mean<0),max_value,network.metabolites(internal),network.actions,[1 .5 .5]); hold on;
im_circle( exp(R_log_mean) .*(R_sign_mean<0),max_value,network.metabolites(internal),network.actions,[1 0 0]);  hold off;
title('Distribution of response coefficients');


% ------------------------------------------------------
% distribution of covariance matrices

C_log_list = log(abs(C_list));
C_log_mean = squeeze(mean(C_log_list,1));
C_log_std  = reshape( sqrt(mean( (reshape(C_log_list,nit,n_int^2) - repmat(reshape(C_log_mean,1,n_int^2),nit,1)).^2 )),n_int,n_int);

C_sign_mean = squeeze(mean(sign(C_list),1));

max_value = max(max(exp(C_log_mean+C_log_std)));

figure(fignum(2)); clf
im_circle( exp(C_log_mean+C_log_std) .*(C_sign_mean>0),max_value,network.metabolites(internal),network.metabolites(internal),[.5 .5 1]); hold on;
im_circle( exp(C_log_mean) .*(C_sign_mean>0),max_value,network.metabolites(internal),network.metabolites(internal),[0 0 1]);  hold on;
im_circle( exp(C_log_mean+C_log_std) .*(C_sign_mean<0),max_value,network.metabolites(internal),network.metabolites(internal),[1 .5 .5]); hold on;
im_circle( exp(C_log_mean) .*(C_sign_mean<0),max_value,network.metabolites(internal),network.metabolites(internal),[1 0 0]);  hold off;
title('Distribution of covariances');


% ------------------------------------------------------
% distribution of correlation matrices

corr_mean = squeeze(mean(corr_list,1));
corr_std  = reshape( sqrt(mean( (reshape(corr_list,nit,n_int^2) - repmat(reshape(corr_mean,1,n_int^2),nit,1)).^2 )),n_int,n_int);

corr_sign_mean = squeeze(mean(sign(corr_list),1));

max_value = 1;

figure(fignum(3)); clf
im_circle( (abs(corr_mean)+corr_std) .*(corr_sign_mean>0),max_value,network.metabolites(internal),network.metabolites(internal),[.5 .5 1]); hold on;
im_circle( abs(corr_mean) .*(corr_sign_mean>0),max_value,network.metabolites(internal),network.metabolites(internal),[0 0 1]);  hold on;
im_circle( (abs(corr_mean)+corr_std) .*(corr_sign_mean<0),max_value,network.metabolites(internal),network.metabolites(internal),[1 .5 .5]); hold on;
im_circle( abs(corr_mean) .*(corr_sign_mean<0),max_value,network.metabolites(internal),network.metabolites(internal),[1 0 0]);  hold off;
title('Distribution of correlations');
