function result = stm_model_to_sampling(network, network_CoHid, v, target_reaction, ind_ignore, print_graphics, psfile_dir, graphics_file)

%STM_MODEL_TO_SAMPLING - Wrapper function for stm_reference_state plus some graphics routines
%
% result = stm_model_to_sampling(network, network_CoHid, v, target_reaction, ind_ignore, print_graphics, psfile_dir, graphics_file)
%
% Calculate a reference metabolic state and show graphics
%
% This function is a wrapper for stm_reference_state and some following graphics routines
%
% Function arguments:
%  network         network data structure
%  network_CoHid   network data structure for graphics
%  v               flux distribution (must be thermodynamically feasible)
%  target_reaction name of target reaction (for synergy calculation)
%  ind_ignore      indices of reactions to be ignored
%  print_graphics  flag
%  psfile_dir      directory for graphics files (only needed if print_graphics = 1)
%  graphics_file   base name for graphics files (only needed if print_graphics = 1)

eval(default('psfile_dir','''''','graphics_file',''''''));

[es_options, es_constraints] = es_default_options(network);

es_options.ind_ignore      = ind_ignore;
es_options.sampling_method = 'accept_flux';
es_options.print_graphics  = print_graphics;
es_constraints.v_mean      = v;

filenames.psfile_dir  = psfile_dir;
filenames.psfile_base = graphics_file;

if ~print_graphics, psfile_dir = []; end

[result, fluxes, es_options, es_constraints] = stm_reference_state(network, es_options, es_constraints);

stm_reference_state_graphics(network, es_options, result, target_reaction, filenames, network_CoHid, [], psfile_dir);

stm_synergy_graphics(network, network_CoHid, result, filenames, target_reaction, psfile_dir);
