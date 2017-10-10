function result = es_model_to_sampling(network, network_CoHid, v, target_reaction, ind_ignore, print_graphics, psfile_dir, graphics_file)

%ES_MODEL_TO_SAMPLING - Wrapper function for es_reference_state plus some graphics routines
%
% result = es_model_to_sampling(network, network_CoHid, v, target_reaction, ind_ignore, print_graphics, psfile_dir, graphics_file)
%
% Calculate a reference metabolic state and show graphics
% 
% This function is a wrapper for es_reference_state and some following graphics routines
%
% Function arguments:
%  network         network data structure
%  network_CoHid   network data structure for graphics
%  v               flux distribution (must be thermodynamically feasible)
%  target_reaction name of target reaction (for synergy calculation)
%  ind_ignore      indices of reactions to be ignored 
%  print_graphics  flag
%  psfile_dir      directory name 
%  graphics_file   file name

eval(default('psfile_dir','''''','graphics_file',''''''));

[es_options, es_constraints] = es_default_options(network);

es_options.ind_ignore      = ind_ignore;
es_options.sampling_method = 'accept_flux';
es_options.print_graphics  = print_graphics;
es_constraints.v_mean      = v;

filenames.psfile_dir           = psfile_dir; 
filenames.reference_state_file = graphics_file;

[result, fluxes, es_options, es_constraints] = es_reference_state(network, es_options, es_constraints);

sampling_reference_state_graphics;

interaction_analysis_graphics(network, network_CoHid, result, filenames, filenames.reference_state_file, target_reaction, es_options.print_graphics,filenames.psfile_dir);
