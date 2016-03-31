function [Sigma_c, Sigma_j, Sigma_c_list, Sigma_j_list, Sigma_c_specific_1Hertz, Sigma_j_specific_1Hertz] = energetics_control_and_noise(network,network_CoHid, c, v, mu, omega_list, tau_list, omega_fast, omega_slow, psfile_dir, basename, es_options, es_constraints, volume)

%[Sigma_c, Sigma_j, Sigma_c_list, Sigma_j_list, Sigma_c_specific_1Hertz, Sigma_j_specific_1Hertz] = energetics_control_and_noise(network,network_CoHid, v, mu, omega_list, tau_list, psfile_dir, basename, es_options, es_constraints, volume)
%
% Analyse control coefficients and noise levels (assuming uniform white noise in each reaction!)
%
% Recommended units: 
%  Enzyme levels in  mol
%  Volume            m^3
%  Metabolite levels mM
%  Fluxes            mol/s
%  Accordingly, the spectral densities will be in mM^2 and (mol/s)^2
% 
% Spectral densities are computed under the assumption that white noise 
% has a spectral density of 1. If white noise is defined as having a covariance 
% given by the Dirac delta function (with prefactor 1), as assumed in the 
% chemical Langevin equation, this corresponds to a Fourier transformation 
% with prefactor 1 for the transformation omega -> t and prefactor 1/2pi for the 
% transformation t -> omega

eval(default('es_options','[]', 'es_constraints','[]','volume','1'));

[nm,nr] = size(network.N);

[K, L, N0, G, pinv_N_R, indep] = network_analyse(network);


% prepare elasticity sampling

if isempty(es_options),         es_options = struct; end
if isempty(es_constraints), es_constraints = struct; end

[es_def_options, es_def_constraints] = es_default_options(network);

% default: simple choice -> keq = 1
es_def_options.seed              = nan;
es_def_constraints.log_c_fix     = 1/RT * es_def_constraints.mu_fix; 
es_def_options.set_alpha_to_half = 1; % set enzymes to half-saturation

es_options            = join_struct(es_def_options,es_options);
es_constraints        = join_struct(es_def_constraints,es_constraints);
es_constraints.mu_fix = mu;
es_constraints.v_fix  = v;
es_options.flag_second_order = 0;


% run elasticity sampling

result = es_sample_model(network.N,network.regulation_matrix,find(network.external),es_constraints,es_options);

% v_molnumbers contains the reaction rates measured in molecules/s
% whereas previous rates must be amounts/time, i.e. mol/s

N_avogadro = 6 * 10^23;

display(sprintf('Assumed volume: %f um^3', 10^18 * volume));
v_molnumbers_plus  = result.v_plus  * N_avogadro;
v_molnumbers_minus = result.v_minus * N_avogadro;
%% conc: values assuming that fluxes are expressed in mM/s
Ec  = result.elasticities.un_E_c;
Ep  = 1/N_avogadro * [diag(sqrt(v_molnumbers_plus)), diag(sqrt(v_molnumbers_plus))];
noise_production   = [v_molnumbers_plus + v_molnumbers_minus] / N_avogadro;

% static contol coefficients

[CJ, CS] = control_coefficients(1/volume*network.N, Ec, network.external);

CJ_sc = diag(1./result.v) * CJ * diag(result.v);
CS_sc = diag(1./result.c) * CS * diag(result.v);


% origin of chemical noise (spectral density of initial rate fluctuations)

% propagation of chemical noise

Sigma_c = [];
Sigma_j = [];
Sigma_c_list = {};
Sigma_j_list = {};
z       = 1;
ind_int = find(network.external ==0);

% The spectral power density is defined as the Fourier transform of the covariance
% function, with the convention 

for omega = omega_list,
  [RS, RJ]     = spectral_response_coefficients(1/volume*N0, L, [], struct('omega',omega),Ec(:,ind_int),Ep); 
  % here: count flux changes relative to flux direction
  RJ = diag(sign(v)) * RJ;
  Sigma_c      = [Sigma_c,      real(diag(RS * RS'))];
  Sigma_j      = [Sigma_j,      real(diag(RJ * RJ'))];
  Sigma_c_list = [Sigma_c_list,{real(RS * RS')}];
  Sigma_j_list = [Sigma_j_list,{real(RJ * RJ')}];
  z = z+1;
end


%--------------------------------------------------------
% GRAPHICS

ca


% -------------------------------------------------------
% static control coefficients

figure(1001); clf; set(gca,'Fontsize',18); 
im(CJ_sc,1.01*max(abs(CJ_sc(:)))); colormap([0.7 0.7 0.7; rb_colors]); 
xlabel('Reaction perturbed'); ylabel('Reaction'); colorbar

figure(1002); clf; set(gca,'Fontsize',18); 
im(CS_sc,1.01*max(abs(CS_sc(:)))); colormap([0.7 0.7 0.7; rb_colors]);
xlabel('Reaction perturbed'); ylabel('Metabolite'); colorbar


% -------------------------------------------------------
% fluctuations caused by chemical noise at omega = 2 pi / s 
% (oscillations of 1 Hertz)

[RS, RJ] = spectral_response_coefficients(N0, L, [], struct('omega',2*pi),Ec(:,ind_int),Ep); 
% here: count flux changes relative to flux direction
RJ = diag(sign(v)) * RJ;

Sigma_c_specific_1Hertz = abs(RS(:,1:nr).^2) + abs(RS(:,nr+1:end).^2);
Sigma_j_specific_1Hertz = abs(RJ(:,1:nr).^2) + abs(RJ(:,nr+1:end).^2);

figure(2001); clf; 
subplot('Position',[0.2 0.1 0.77 0.7]);%set(gca,'Fontsize',18); 
im(Sigma_j_specific_1Hertz,[],network.actions,network.actions); colormap(rb_colors); 
axis square; my_xticklabel; xlabel('Reaction (origin of chemical noise)'); colorbar

figure(2002); clf; 
subplot('Position',[0.4 0.1 0.57 0.7]);
im(Sigma_c_specific_1Hertz,[],network.metabolites(ind_int),network.actions); colormap(rb_colors); 
my_xticklabel; xlabel('Reaction (origin of chemical noise)'); colorbar


% -------------------------------------------------------
% display spectral densities

figure(2); clf; 
hold on; 
h = plot(1/(2*pi) * omega_list,Sigma_c'); 
set(gca,'XScale','Log', 'YScale','Log');set(gca,'Fontsize',26) 
axis tight; line_colors(h,'copper'); 
%legend(escape_underscores(network.metabolites(ind_int)),'Fontsize',6);
xlabel('Frequency [1/s]'); ylabel('Spectral density [mM^2 s]');

figure(3); clf; 
hold on; 
h = plot(1/(2*pi) * omega_list,Sigma_j'); set(gca,'XScale','Log', 'YScale','Log'); 
axis tight; line_colors(h,'copper'); set(gca,'Fontsize',26);
xlabel('Frequency [1/s]'); ylabel('Spectral density [(mol/s)^2 s]');
%legend(escape_underscores(network.actions),'Fontsize',6,'Location','EastOutside');


%figure(4); plot(es_constraints.log_c_fix,'-*');

dum = nan*ones(nm,1);

% -------------------------------------------------
% Spectral densities

%figure(10); dum(ind_int) =  Sigma_c(:,1)'; netgraph_molnumbersentrations(network_CoHid, dum,[],0,struct('colormap',rb_colors)); title('Concentrations, slow fluctuations');

%figure(11); dum(ind_int) =  Sigma_c(:,end)'; netgraph_molnumbersentrations(network_CoHid, dum,[],0,struct('colormap',rb_colors)); title('Concentrations, fast fluctuations');

%figure(12); netgraph_molnumbersentrations(network_CoHid, [], Sigma_j(:,1)',0,struct('arrowstyle','none','colormap',rb_colors));   title('Fluxes, slow fluctuations');

%figure(13); netgraph_molnumbersentrations(network_CoHid, [], Sigma_j(:,end)',0,struct('arrowstyle','none','colormap',rb_colors)); title('Fluxes, fast fluctuations');


% -------------------------------------------------
% Spectral metabolite correlations (threshold for display: 10%)

dum = zeros(nm,nm);
dum(ind_int,ind_int) = Sigma_c_list{ find(omega_list == omega_slow) }; 
figure(110); 
netgraph_metabolite_interactions(network_CoHid, diag(dum), dum, rb_colors, struct('relative_threshold',0.1)); 

% 1/sec
dum(ind_int,ind_int) = Sigma_c_list{ find(omega_list == omega_fast) };  
figure(111);
netgraph_metabolite_interactions(network_CoHid, diag(dum), dum, rb_colors, struct('relative_threshold',0.1)); 


% -------------------------------------------------
% Spectral flux correlations (threshold for display: 10%)

gp = struct('relative_threshold',0.1,'arrowstyle','none','actprintnames',1,'metprintnames',0);

% 1/(1000) 1/sec  (16.7 minutes period)
dum = Sigma_j_list{ find(omega_list == 2 * pi * 1/(1000)) };
figure(112); interaction_network_plot(network_CoHid, diag(dum), dum, rb_colors, gp); 

% 1/sec
dum = Sigma_j_list{ find(omega_list == 2 * pi * 1) };
figure(113); interaction_network_plot(network_CoHid, diag(dum), dum, rb_colors, gp); 

% -------------------------------------------------

figure(20); netgraph_concentrations(network_CoHid, [],abs(result.A),1,struct('colormap',rb_colors,'arrowstyle','none')); 

figure(21); netgraph_concentrations(network_CoHid, result.mu,[],1,struct('colormap',rb_colors)); 

figure(22); netgraph_concentrations(network_CoHid, [] , 1/max(noise_production) * noise_production,1,struct('colormap',rb_colors)); 

% -------------------------------------------------
% variance quantities measured on different time scales

tau_list = 10.^[-2:0.1:4];
Var_c = spectral_density_to_smooth_variance(omega_list,Sigma_c,tau_list);
Var_j = spectral_density_to_smooth_variance(omega_list,Sigma_j,tau_list);

figure(1011); clf; 
h = plot(tau_list,sqrt(Var_c)); 
set(gca,'XScale','log','YScale','log','Fontsize',26);
xlabel('Observation time scale \tau [s]'); ylabel('Concentration Std Dev [mM]'); 
line_colors(h,'copper')
axis tight
%legend(escape_underscores(network.metabolites(network.external==0)),'Fontsize',6);

figure(1012); clf; 
h = plot(tau_list,sqrt(Var_j)); 
set(gca,'XScale','log','YScale','log','Fontsize',26)
xlabel('Observation time scale \tau [s]'); ylabel('Flux Std Dev [mol/s]');
line_colors(h,'copper')
axis tight
%legend(escape_underscores(network.actions),'Fontsize',6);

figure(1013); clf; set(gca,'Fontsize',18)
h = plot(tau_list,diag(1./c(find(network.external==0))) * sqrt(Var_c)); set(gca,'XScale','log','YScale','log');
xlabel('Observation time scale \tau [s]'); ylabel('Concentration CV [unitless]'); 
line_colors(h,'copper')
%legend(escape_underscores(network.metabolites(network.external==0)),'Fontsize',6);

figure(1014); clf; set(gca,'Fontsize',18)
h = plot(tau_list,diag(1./v)*sqrt(Var_j)); set(gca,'XScale','log','YScale','log');
xlabel('Observation time scale \tau [s]'); ylabel('Flux CV [unitless]');
line_colors(h,'copper')
%legend(escape_underscores(network.actions),'Fontsize',6);


if length(psfile_dir),
cd(psfile_dir);
display(sprintf('Saving graphics to directory %s', psfile_dir));
print([ psfile_dir '/' basename '_flux_control_coeff.eps'], '-f1001', '-depsc');
print([ psfile_dir '/' basename '_conc_control_coeff.eps'], '-f1002', '-depsc');
print([ psfile_dir '/' basename '_flux_chem_varation.eps'], '-f2001', '-depsc');
print([ psfile_dir '/' basename '_conc_chem_varation.eps'], '-f2002', '-depsc');
print([ psfile_dir '/' basename '_spectral_conc_freq.eps'], '-f2', '-depsc');
print([ psfile_dir '/' basename '_spectral_flux_freq.eps'], '-f3', '-depsc');
%print([ psfile_dir '/' basename '_spectral_conc_slow_network.eps'], '-f10', '-depsc');
%print([ psfile_dir '/' basename '_spectral_conc_fast_network.eps'], '-f11', '-depsc');
%print([ psfile_dir '/' basename '_spectral_flux_slow_network.eps'], '-f12', '-depsc');
%print([ psfile_dir '/' basename '_spectral_flux_fast_network.eps'], '-f13', '-depsc');
print([ psfile_dir '/' basename '_spectral_corr_conc_slow_network.eps'], '-f110', '-depsc');
print([ psfile_dir '/' basename '_spectral_corr_conc_fast_network.eps'], '-f111', '-depsc');
print([ psfile_dir '/' basename '_spectral_corr_flux_slow_network.eps'], '-f112', '-depsc');
print([ psfile_dir '/' basename '_spectral_corr_flux_fast_network.eps'], '-f113', '-depsc');
print([ psfile_dir '/' basename '_spectral_affinities_network.eps'], '-f20', '-depsc');
print([ psfile_dir '/' basename '_spectral_chem_pot_network.eps'], '-f21', '-depsc');
print([ psfile_dir '/' basename '_spectral_noise_production_network.eps'], '-f22', '-depsc');
print([ psfile_dir '/' basename '_spectral_conc_smoothed_variance.eps'], '-f1011', '-depsc');
print([ psfile_dir '/' basename '_spectral_flux_smoothed_variance.eps'], '-f1012', '-depsc');
end

figure(1001); title('Scaled flux control coefficients');          
figure(1002); title('Scaled concentration control coefficients'); 
figure(2001); xlabel('Flux spectral density [(mol/s)^2 s] at 1/s');          
figure(2002); xlabel('Concentration spectral density [mM^2 s] at 1/s');          
figure(2);    title('Spectral densities of concentrations');
figure(3);    title('Spectral densities of reaction rates');
figure(110);  title('Concentrations, slow fluctuations');
figure(111);  title('Concentrations, fast fluctuations');
figure(112);  title('Fluxes, slow fluctuations');
figure(113);  title('Fluxes, fast fluctuations');
figure(20);   title('Absolute reaction affinities');
figure(21);   title('Chemical potentials');
figure(22);   title('Noise production');
