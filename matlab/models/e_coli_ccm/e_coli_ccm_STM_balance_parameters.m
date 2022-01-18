% -------------------------------------------------
% Determine balanced parameter set for E coli model and save to file
% /resource-data/escherichia_coli/data-kinetic/e_coli_balanced_parameters.tsv
% -------------------------------------------------

clear; 

save_files = 0;

model_name              = 'e_coli_noor_2016_glucose';
DATA_DIR                = [es_BASEDIR  '/../resource-data/models-article/Escherichia_coli_ccm/original_model/'];
kinetic_data_file       = [DATA_DIR 'e_coli_noor_2016_kinetic_data_modified.tsv'];

% deprecated:
%model_file              = [ DATA_DIR 'e_coli_noor_2016_plus_overflow_oxa_internal.tsv'];
%balanced_parameter_file = [ DATA_DIR 'e_coli_noor_2016_plus_overflow_oxa_internal_balanced_kinetics.tsv'];

model_file              = [ DATA_DIR 'other/e_coli_noor_2016_plus_overflow.tsv'];
balanced_parameter_file = [ DATA_DIR 'other/e_coli_noor_2016_plus_overflow_balanced_kinetics.tsv'];

prior_file = [es_BASEDIR  '../resource-data/data/data-general/cmb_prior.tsv'];

pb_options                = parameter_balancing_options;
pb_options.use_sbml_ids   = 0;
pb_options.use_kegg_ids   = 1;
balanced_parameters_SBtab = parameter_balancing(model_file, balanced_parameter_file, kinetic_data_file, prior_file, [], model_name, pb_options);

if save_files, 
  print([DATA_DIR '/e_coli_balanced_parameters_fig_1.eps'],'-f1','-depsc');
  print([DATA_DIR '/e_coli_balanced_parameters_fig_2.eps'],'-f2','-depsc');
end
