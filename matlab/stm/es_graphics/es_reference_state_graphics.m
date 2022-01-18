function es_reference_state_graphics(network, es_options, result, target_reaction, es_filenames, network_CoHid, network_CoSplit, psfile_dir)

% -------------------------------------------------------------
% graphics for elasticity sampling reference states

eval(default('network_CoHid','network','network_CoSplit','[]'));
if isempty(network_CoSplit),
  network_CoSplit = network_CoHid;
end

eval(default('target_reaction','''Biomass production''','target_product','''Biomass''')); % 'Maintenance'
if exist('es_filenames','var'),
  eval(default('psfile_dir','es_filenames.psfile_dir'));
else
  eval(default('psfile_dir','[]'));
end
  
N       = network.N;
[nm,nr] = size(N);
ind_ext = find(network.external);

p.print = es_options.print_graphics;  % use layout for printing

% network
figure(1); es_graphics(network,network_CoSplit,'names',struct('actprintnames',0,'FontSize',24));

% figure(1000); clf; 
%es_graphics(network,network_CoHid,'fluxes',struct('v',fluxes.v_mean),p);

%figure(2); clf; es_graphics(network,network_CoHid,'fluxes',struct('v',result.v,'production_rates',N*result.v),join_struct(p,struct('show_regulation',0)));
figure(2); clf; es_graphics(network,network_CoHid,'fluxes',struct('v',result.v),join_struct(p,struct('show_regulation',0,'arrowsize',0.05)));
%figure(2); title(sprintf('Fluxes and net production rates')); 

% flux movies
% figure(102); clf; M = es_graphics(network,network_CoHid,'flux_movie',struct('v',fluxes.best_v,'n_frames',6),p);
% movie(M,1,4);

figure(3); clf; es_graphics(network,network_CoSplit,'chemical_potentials',struct('mu',result.mu),p);

figure(6); clf; es_graphics(network,network_CoHid,'thermodynamic_forces',struct('theta',result.A/RT,'v',result.v),p);

figure(7); clf; es_graphics(network,network_CoSplit,'elasticities',struct('elasticities',result.saturation.beta_M + result.saturation.beta_I + result.saturation.beta_A,'v',result.v),p);

figure(8); clf; es_graphics(network,network_CoSplit,'elasticities',struct('elasticities',result.elasticities.sc_E_c,'v',result.v),p);

figure(9); clf; es_graphics(network,network_CoHid,'dissipation',struct('theta',result.A/RT,'v',result.v),p);

figure(10); clf; es_graphics(network,network_CoSplit,'metabolites',struct('c',result.c),p);

figure(11); clf; es_graphics(network,network_CoHid,'enzymes',struct('u',result.u .* result.KV),p);

if length(target_reaction),
  display(sprintf('Showing response and synergies w.r.t. target reaction %s',target_reaction));
  %% response coefficients for target reaction
  ind_target = find(strcmp(network.actions,target_reaction));
  if isempty(ind_target), error(sprintf('Unknown target reaction %s', target_reaction)); end
  rs_enz     = column(result.control.RJu_sc(ind_target,:));
  rs_ext     = zeros(nm,1); 
  rs_ext(ind_ext) = result.control.RJs_sc(ind_target,:);
  output = network.actions{ind_target};
else,
  %% response coefficients production of (external) target metabolite
  display(sprintf('Showing response and synergies w.r.t. target metabolite %s',target_product));
  ind_target    = find(strcmp(network.metabolites,target_product));
  if network.external(ind_target)==0, error('Production of internal metabolite makes no sense'); end 
  rs_enz = column(N(ind_target,:) * result.control.RJu_sc);
  rs_ext = zeros(nm,1); rs_ext(ind_ext) = N(ind_target,:) * result.control.RJs_sc;
  output = [network.metabolites{ind_target} ' production']; 
end

rs_enz(ind_target,1) = nan;
figure(4); clf; 
es_graphics(network,network_CoHid,'response_coefficients',struct('R',rs_enz,'v',result.J),p);

figure(5); clf; 
es_graphics(network,network_CoSplit,'response_coefficients_Sext',struct('R',rs_ext,'v',result.J),p);


% -------------------------------------------------------------
% save graphics

if es_options.print_graphics,
  cd(psfile_dir);
  display(sprintf('Saving graphics to directory %s',psfile_dir));
  print( [ es_filenames.psfile_base '_network.eps'],   '-f1', '-depsc');
  print( [ es_filenames.psfile_base '_fluxes.eps'],    '-f2', '-depsc');
  print( [ es_filenames.psfile_base '_chempot.eps'],   '-f3', '-depsc');
  print( [ es_filenames.psfile_base '_resp_enz.eps'],  '-f4', '-depsc');
  print( [ es_filenames.psfile_base '_resp_ext.eps'],  '-f5', '-depsc');
  print( [ es_filenames.psfile_base '_forces.eps'], '-f6', '-depsc');
  print( [ es_filenames.psfile_base '_saturation.eps'], '-f7', '-depsc');
  print( [ es_filenames.psfile_base '_elasticities.eps'],'-f8', '-depsc');
  print( [ es_filenames.psfile_base '_dissipation.eps'],'-f9', '-depsc');
  print( [ es_filenames.psfile_base '_concentrations.eps'],'-f10', '-depsc');
  print( [ es_filenames.psfile_base '_vmax_geom.eps'],'-f11', '-depsc');
  %%  movie_save([ result_file '_fluxes'],M);
end

% figure(1); title(sprintf('Measured fluxes and external metabolites')); 
figure(2); title(sprintf('Fluxes and net production rates')); 
figure(3); title(sprintf('Chemical potentials'));
figure(4); title(sprintf('Scaled resp coeff. (on output %s\n by enz)',escape_underscores(output))); 
figure(5); title(sprintf('Scaled resp coeff. (on output %s\n ext metab)',escape_underscores(output)));
figure(6); title(sprintf('Driving forces (in kJ/mol)'));
figure(7); title(sprintf('Saturation values'));
figure(8); title(sprintf('Scaled elasticities'));
figure(9); title(sprintf('GFE dissipation'));
figure(10); title(sprintf('log10 concentrations'));
figure(11); title(sprintf('log10 Vmax geometric means')); drawnow


% -----------------------------------------------------------
% distribution of elasticities

% Cumulative distributions: 

% distribution of elasticities (sorted values)
% figure(100); plot(sort(result.elasticities.sc_E_c(result.elasticities.sc_E_c~=0)))
% title('Distribution of scaled elasticities (sorted values)');

% figure(101); plot(sort(result.elasticities.un_E_c(result.elasticities.un_E_c~=0)))
% title('Distribution of unscaled elasticities (sorted values)');

figure(202);  clf
plot(result.elasticities.sc_E_c(:),result.elasticities.un_E_c(:),'o');
axis tight;
xlabel('Scaled elasticities');
ylabel('Unscaled elasticities');
title('Elasticities')