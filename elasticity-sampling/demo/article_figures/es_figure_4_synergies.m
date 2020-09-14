% Run elasticity sampling and analyse synergies 
%
% See Article Figure 5


% -------------------------------------------------------------
% load model and choose reference state

model_dir  = [es_BASEDIR '/../resources/models-article/Example_Loop_Pathway'];
model_name = 'Example_Loop_Pathway'; 
ref_state  = 'Example_Loop_Pathway_reference_state'; 

kinetic_law = 'cs';
nrun        = 50; 
target_flux = 'E6';

cd(model_dir); load(model_name); load(ref_state);

c = reference_state.c;
v = reference_state.v; 

dmu = -reference_state.A;    
mu  = [];    
v   = 1/max(abs(v)) * v;


% -------------------------------------------------------------
% define network structure: matrices N and W, ind_ext, ind_target_flux

if ~exist('network_CoHid','var'), network_CoHid = network; end

N               = network.N;
[nm,nr]         = size(N);
W               = network.regulation_matrix;
ind_ext         = find(network.external);
ind_target_flux = find(strcmp(target_flux,network.actions));


% -------------------------------------------------------------
% settings for graphics 

gpp = struct('flag_edges',1, 'squaresize',0.1,'arrowsize',0.05,'actprintnames',1,'metstyle','fixed','metvalues',zeros(size(network.metabolites)),'actstyle','none','arrowvalues',v, 'linecolor',[0 0 0],'arrowcolor',[.7 .7 .7],'colorbar',0, 'FontSize',16,'text_offset',[.04,.05],'colormap',rb_colors,'relative_threshold',0.01);


%% note that values in interaction graphics are thresholded!

gp = gpp;

figure(1002); clf; 
netgraph_concentrations(network_CoHid, c, v, 1, gp)


% -------------------------------------------------------------
% define constraints for sampling (concentrations in mM, fluxes in mM/s, energies in kJ/mol)

if ~exist('es_options'),
  [es_options,es_constraints] = es_default_options(network);
end

es_options.kinetic_law     = kinetic_law;
es_options.sampling_method = 'accept flux';
es_constraints.v_fix       = v;
if length(dmu), es_constraints.dmu_fix = dmu; end

% -------------------------------------------------------------
% Single sampling run

es_options.seed              = 0;
es_options.set_alpha_to_half = 1;

reference_state = es_sample_model(N,W,ind_ext,es_constraints,es_options);
RFuu   = squeeze(reference_state.control.RJuu_sc(ind_target_flux,:,:));

RFu = reference_state.control.RJu_un(ind_target_flux,:);


% ------------------------------------
% For comparison, run fba (all fluxes between -1 and 1), optimise benefit target flux

fba_constraints                     = fba_default_options(network);
fba_constraints.zv                  = zeros(size(network.actions));
fba_constraints.zv(ind_target_flux) = 1;

inhibition_factor = 0.9;

[fba_a_f_ref,fba_a_f_single,fba_a_f_double]    = fba_interactions(network,fba_constraints,inhibition_factor);
[moma_a_f_ref,moma_a_f_single,moma_a_f_double] = moma_interactions(network,fba_constraints,inhibition_factor,v);

epistasis_a            = log_epistasis(fba_a_f_ref,fba_a_f_single,fba_a_f_double);
epistasis_moma_a       = log_epistasis(moma_a_f_ref,moma_a_f_single,moma_a_f_double);
segre_epistasis_a      = segre_epistasis(fba_a_f_ref,fba_a_f_single,fba_a_f_double);
segre_epistasis_moma_a = segre_epistasis(moma_a_f_ref,moma_a_f_single,moma_a_f_double);

inhibition_factor = 0.5;

[fba_b_f_ref,fba_b_f_single,fba_b_f_double]    = fba_interactions(network,fba_constraints,inhibition_factor);
[moma_b_f_ref,moma_b_f_single,moma_b_f_double] = moma_interactions(network,fba_constraints,inhibition_factor,v);

epistasis_b            = log_epistasis(fba_b_f_ref,fba_b_f_single,fba_b_f_double);
epistasis_moma_b       = log_epistasis(moma_b_f_ref,moma_b_f_single,moma_b_f_double);
segre_epistasis_b      = segre_epistasis(fba_b_f_ref,fba_b_f_single,fba_b_f_double);
segre_epistasis_moma_b = segre_epistasis(moma_b_f_ref,moma_b_f_single,moma_b_f_double);

inhibition_factor = 0.1;

[fba_c_f_ref,fba_c_f_single,fba_c_f_double]    = fba_interactions(network,fba_constraints,inhibition_factor);
[moma_c_f_ref,moma_c_f_single,moma_c_f_double] = moma_interactions(network,fba_constraints,inhibition_factor,v);

epistasis_c            = log_epistasis(fba_c_f_ref,fba_c_f_single,fba_c_f_double);
epistasis_moma_c       = log_epistasis(moma_c_f_ref,moma_c_f_single,moma_c_f_double);
segre_epistasis_c      = segre_epistasis(fba_c_f_ref,fba_c_f_single,fba_c_f_double);
segre_epistasis_moma_c = segre_epistasis(moma_c_f_ref,moma_c_f_single,moma_c_f_double);


% ------------------------------------
% steady state concentrations and fluxes 

gp = gpp;

figure(1); clf; 
netgraph_concentrations(network_CoHid,reference_state.c,reference_state.J,1,gp);

% plot colorbar
figure(1001); 
clf; my_colorbar(0,max(reference_state.c),5,'east',0,rb_colors); axis off; 


% ------------------------------------
% second flux response coefficients (matrix)

inhibition_factor = 0.1;

RFuu_sc = RFuu * [1-inhibition_factor]^2;

figure(2); clf; set(gca,'FontSize',20);
mmax = max(max(abs(triu(RFuu_sc,1))));
im(RFuu_sc,[mmax],network.actions);%,network.actions); 
colormap(rb_colors); axis equal; axis tight; colorbar


% ------------------------------------
% second flux response coefficients (on network)

gp = gpp;
gp.relative_threshold = 0.01;

figure(3); clf; 
interaction_network_plot(network_CoHid,[], sign(RFuu), rb_colors, gp);

figure(4); clf;
interaction_network_plot(network_CoHid,[], RFuu, rb_colors, gp)


% -------------------------------------------------------------
% show fba results 

figure(2001); clf; set(gca,'FontSize',20);
im(epistasis_a,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(2002); clf;  set(gca,'FontSize',20);
im(epistasis_b,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(2003); clf;  set(gca,'FontSize',20);
im(epistasis_c,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(2011); clf;  set(gca,'FontSize',20);
im(epistasis_moma_a,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(2012); clf;  set(gca,'FontSize',20);
im(epistasis_moma_b,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(2013); clf;  set(gca,'FontSize',20);
im(epistasis_moma_c,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(2101); clf; 
interaction_network_plot(network_CoHid, [], epistasis_a, rb_colors, gpp);

figure(2102); clf;
interaction_network_plot(network_CoHid, [], epistasis_b, rb_colors, gpp);

figure(2103); clf;
interaction_network_plot(network_CoHid, [], epistasis_c, rb_colors, gpp);

figure(2111); clf;
interaction_network_plot(network_CoHid, [], epistasis_moma_a, rb_colors, gpp);

figure(2112); clf;
interaction_network_plot(network_CoHid, [], epistasis_moma_b, rb_colors, gpp);

figure(2113); clf;
interaction_network_plot(network_CoHid, [], epistasis_moma_c, rb_colors, gpp);

figure(3001); clf; set(gca,'FontSize',20);
im(segre_epistasis_a,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(3002); clf;  set(gca,'FontSize',20);
im(segre_epistasis_b,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(3003); clf;  set(gca,'FontSize',20);
im(segre_epistasis_c,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(3011); clf;  set(gca,'FontSize',20);
im(segre_epistasis_moma_a,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(3012); clf;  set(gca,'FontSize',20);
im(segre_epistasis_moma_b,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(3013); clf;  set(gca,'FontSize',20);
im(segre_epistasis_moma_c,[],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar

figure(3101); clf;
interaction_network_plot(network_CoHid, [], segre_epistasis_a, rb_colors, gpp);

figure(3102); clf; 
interaction_network_plot(network_CoHid, [], segre_epistasis_b, rb_colors, gpp);

figure(3103); clf; 
interaction_network_plot(network_CoHid, [], segre_epistasis_c, rb_colors, gpp);

figure(3111); clf; 
interaction_network_plot(network_CoHid, [], segre_epistasis_moma_a, rb_colors, gpp);

figure(3112); clf; 
interaction_network_plot(network_CoHid, [], segre_epistasis_moma_b, rb_colors, gpp);

figure(3113); clf; 
interaction_network_plot(network_CoHid, [], segre_epistasis_moma_c, rb_colors, gpp);

% -------------------------------------------------------------
% many runs

es_options.set_alpha_to_half = 0;

for it = 1:nrun,
  display(sprintf('%d/%d',it,nrun));
  es_options.seed        = it;
  reference_statelist{it}      = es_sample_model(N,W,ind_ext,es_constraints,es_options);
  clist(:,it)         = reference_statelist{it}.c;
  vlist(:,it)         = reference_statelist{it}.v;
  RFu_sc_list(it,:)    = reference_statelist{it}.control.RJu_sc(ind_target_flux,:);
  RFuu_sc_list(it,:,:) = reference_statelist{it}.control.RJuu_sc(ind_target_flux,:,:);
end

c_mean      = c;
v_mean      = v;
RFuu_sc_frac = squeeze(mean(sign(RFuu_sc_list),1));
RFuu_sc_mean = squeeze(mean(RFuu_sc_list,1));
RFuu_sc_std  = squeeze(std(RFuu_sc_list));
RFu_sc_mean  = mean(RFu_sc_list,1);
RFu_sc_std   = std(RFu_sc_list);


% ---------------------------------------------------
% second flux response coefficients

% second flux response coefficients; fraction of positive/negative signs

figure(12); clf;  
set(gca,'FontSize',20);
im(RFuu_sc_frac,[-1,1],network.actions);
colormap(rb_colors); axis equal; axis tight

gp = gpp;
gp.relative_threshold = 0.01;

figure(13); clf; 
interaction_network_plot(network_CoHid, [], RFuu_sc_frac, rb_colors, gp); 

figure(1013); clf; 
my_colorbar(-1,1,4,'east',0,rb_colors); axis off; 

% second flux response coefficients; mean values

inhibition_factor = 0.9;
RFuu_sc_mean_sc = RFuu_sc_mean * [1-inhibition_factor]^2;

figure(15); clf;
interaction_network_plot(network_CoHid, [], RFuu_sc_mean_sc, rb_colors, gp); 


% -------------------------------------------------------------
% second flux response coefficients

figure(14); clf; set(gca,'FontSize',20);
mmax = max(max(abs(triu(RFuu_sc_mean,1))));
im(RFuu_sc_mean,[mmax],network.actions);
colormap(rb_colors); axis equal; axis tight; colorbar


% ------------------------------------------------------------
% add titles

figure(1); title('Steady state values');
figure(2); title(sprintf('Unscaled 2^{nd} resp. coeff.: flux[%d]/enzyme ',ind_target_flux));
figure(3); title('Scaled 2nd order resp. coefficients - Signs');
figure(4); title('Scaled 2nd order resp. coefficients - Normalised values');
figure(2001); title('Synergies from FBA (inhibition factor = 0.9)');
figure(2002); title('Synergies from FBA (inhibition factor = 0.5)');
figure(2102); title('Synergies from FBA (inhibition factor = 0.1)');
figure(2003); title('Synergies from FBA (inhibition factor = 0.9)');
figure(2004); title('Synergies from FBA (inhibition factor = 0.5)');
figure(2104); title('Synergies from FBA (inhibition factor = 0.1)');
figure(2005); title('Segre epistasis from FBA analysis (inhibition factor = 0.9)');
figure(2006); title('Segre epistasis from FBA analysis (inhibition factor = 0.5)');
figure(2106); title('Segre epistasis from FBA analysis (inhibition factor = 0.1)');
figure(2007); title('Segre epistasis from FBA analysis (inhibition factor = 0.9)');
figure(2008); title('Segre epistasis from FBA analysis (inhibition factor = 0.5)');
figure(2108); title('Segre epistasis from FBA analysis (inhibition factor = 0.1)');
figure(12); title(sprintf('Sign fraction unscaled 2^{nd} resp. coefficients: flux[%d]/enzyme',ind_target_flux));
figure(13); title('Scaled 2nd order resp. coefficients: fraction(positive)');
figure(15); title('Scaled 2nd order resp. coefficients: Mean values');
figure(14); title(sprintf('Unscaled 2^{nd}  resp. coefficients: flux[%d]/enzyme',ind_target_flux));
