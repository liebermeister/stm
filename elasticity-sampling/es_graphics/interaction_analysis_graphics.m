function interaction_analysis_graphics(network, network_CoHid, result, es_filenames, target_reaction, psfile_dir, R_target_u_sc,R_target_uu_sc,flag_analyse_synergy_cycles,n_cluster);

% interaction_analysis_graphics(network, network_CoHid, result, es_filenames, result_file, target_reaction, print_style, psfile_dir, R_target_u_sc,R_target_uu_sc,flag_analyse_synergy_cycles,n_cluster);
%
% Statistical analysis of first- and second-order response coefficients 
% based a previous single elasticity sampling run
%
% - elasticity sampling of reaction fluxes and chemical potentials
% - statistics over response coefficients
% - compute epistastis measure
%
% read data from [result_file] in directory [es_filenames.data_dir 'es_sampling/']
% save graphics to directory 'psfile_dir'

eval(default('R_target_u_sc','[]','R_target_uu_sc','[]','n_cluster','10','flag_analyse_synergy_cycles','1','psfile_dir','[]'));

if length(psfile_dir),
  network_CoHid.graphics_par.FontSize = 6;
  network.graphics_par.FontSize       = 6;
else
  network_CoHid.graphics_par.FontSize = 10;
  network.graphics_par.FontSize       = 10;
end

% fix (to be removed later)
network_CoHid.graphics_par = join_struct(struct('omitreactions',[],'omitmetabolites',[]),network_CoHid.graphics_par);

% ------------------------------------------------------
% load prepared results

if isempty(network),
  cd([es_filenames.es_dir]); 
  load(result_file);            % load 'options', 'result'
  load(es_filenames.network_file); % load 'network'
end

N       = network.N;
W       = network.regulation_matrix;
[nm,nr] = size(N);
ind_ext = find(network.external);


% -------------------------------------------------------------
% Extract influence (1st order scaled response coefficients) 
% and interaction (2nd order scaled response coefficients)
%  here self-interactions are only set to zero, not to nan
%  make sure they are omitted from statistics later on

if isempty(R_target_u_sc),
  ind_target     = find(strcmp(target_reaction,network.actions));
  R_target_u_sc  = column(result.control.RJu_sc(ind_target,:));
  R_target_uu_sc = squeeze(result.control.RJuu_sc(ind_target,:,:));
  % R_target_u_sc  = result.control.Rtarget_sc_u;
  % R_target_uu_sc = result.control.Rtarget_sc_uu;
end
influence   = R_target_u_sc;
interaction = R_target_uu_sc - diag(diag(R_target_uu_sc));


% --------------------------------------------------------------
% Thresholding (1 percent quantiles) and scaling 
% WITHIN THE SUBNETWORK TO BE SHOWN -> 'interaction_thr'

ind_show = network_CoHid.graphics_par.reaction_mapping;
dum = interaction(ind_show,ind_show); 
dum = dum(find(triu(ones(size(dum,1)))));

q_lo            = quantile(dum,0.01);
q_hi            = quantile(dum,0.99);
interaction_thr = interaction;
interaction_thr(find(double(interaction_thr>q_lo) .* double(interaction_thr<q_hi))) = 0;
interaction_thr = interaction_thr/nanmax(abs(dum)); 
interaction_thr = 0.5*(interaction_thr + interaction_thr');


% --------------------------------------------------------------
% interaction degree: number of interactions per enzyme; 
% 98% quantile as threshold value, leading to a mean degree of about 0.02 * nr

q_thr     = quantile(abs(interaction(:)),0.98);

synergy_degree     = sum(abs(interaction)>q_thr)';
synergy_degree_pos = sum(interaction>q_thr)';
synergy_degree_neg = sum(interaction<-q_thr)';

% use upper five percent for computing sign ratios
q_thr     = quantile(abs(interaction(:)),0.75);
interaction_sign = sign(interaction .* [abs(interaction)>q_thr] );
log2_sign_ratio  = log2(sum([interaction_sign>0]) ./ sum([interaction_sign<0]));

total_positive_synergies = sum(sum(triu(interaction_sign>0)))
total_negative_synergies = sum(sum(triu(interaction_sign<0)))


% -----------------------------------------------------------
% display response coefficients (1st and 2nd order) on network

gp = struct('actstyle','fixed','linecolor',[0 0 0],'arrowcolor',[.7 .7 .7],'colorbar',1,'text_offset',[.01,-.01],'colormap',rb_colors, 'actprintnames',0,'hold_on',1,'FontSize',network.graphics_par.FontSize); 
% 'squaresize',0.03, 'arrowsize',0.03,'arrowstyle', 'fluxes','arrowvalues',result.J, 'arrowvaluesmax',max(abs(result.J)),

figure(101); clf; 
gp.actvaluesmax  = max(abs(influence(network_CoHid.graphics_par.reaction_mapping)));
gp.actprintnames = 0;
netgraph_concentrations(network_CoHid,[], influence',1,gp);

ind_omit = label_names(network_CoHid.graphics_par.omitreactions,network_CoHid.actions);
ind_omit = ind_omit(find(ind_omit));
interaction_thr_show = interaction_thr;
interaction_thr_show(:,ind_show(ind_omit)) = nan; 

figure(102); clf; 
interaction_network_plot(network_CoHid, zeros(size(network.actions)),interaction_thr_show, rb_colors, gp); 

% interaction degrees
figure(103); clf; 
gp.actvaluesmax = max( synergy_degree(network_CoHid.graphics_par.reaction_mapping));
netgraph_concentrations(network_CoHid,[], synergy_degree,1,gp);


% -------------------------------------------------------------
% plot first-order versus second-order response coefficients

CV                 = result.control.RJu_sc;
CV_mean            = 1/2*[CV+CV'];
influence_products = column(influence) * column(influence)';

% replace diagonal elements by nan: 'interactions_nan'
interactions_nan = interaction; 
interactions_nan(find(eye(size(interactions_nan)))) = nan; 

% remove very small interactions
interactions_nan(find(abs(interactions_nan)<10^-3))=nan; 


% ---------------------------------------------------------------
% interactions versus individual influences

ind_triu = find(triu(ones(size(interactions_nan,1)))); 

figure(111); clf; set(gca,'FontSize',18)
hist(interactions_nan(ind_triu), 200)
xlabel('Synergy'); ylabel('Count number (close-up)');
a = axis; axis([a(1) a(2), 0, sqrt(a(4))]);

figure(1111); clf; set(gca,'FontSize',18); 
influences = repmat(influence,length(influence),1);
influences(find(eye(size(influences)))) = nan; % remove diagonal elements
plot(influences(ind_triu),interactions_nan(ind_triu),'.','Color',[.7 .8 1],'MarkerSize',15);  hold on 
axis_almost_tight; a = axis; 
plot([a(1) a(2)],[0 0],'--k'); hold on; 
plot([0 0],[a(3) a(4)],'--k'); hold on;
errorbar(influence,nanmean(interactions_nan),nanstd(interactions_nan),'.','Color',[.2 .2 .8],'LineWidth',2); hold off
xlabel('Influence (1st order scaled control)'); ylabel('Synergy (2nd order scaled control)');
a1 = nanquantile(influences(ind_triu),0.001);
a2 = nanquantile(influences(ind_triu),0.999);
a3 = nanquantile(interactions_nan(ind_triu),0.001);
a4 = nanquantile(interactions_nan(ind_triu),0.999);
if prod(double(isfinite([a1 a2 a3 a4]))),
  axis([a1 a2 a3 a4]);
end

% ----
% interactions versus products of influences

figure(112); clf;
influence_products_nan = influence_products;
influence_products_nan(abs(influence_products_nan)<10^-5) = nan;
set(gca,'FontSize',18)
plot(influence_products(ind_triu),interactions_nan(ind_triu),'.','Color',[.2 .2 0.8],'MarkerSize',15);  hold on 
axis_almost_tight; a = axis;
plot([a(1) a(2)],[0 0],'--k'); hold on; plot([0 0],[a(3) a(4)],'--k'); hold on;
xlabel('Influence product'); ylabel('Synergy');

figure(1112); clf; set(gca,'FontSize',18)
hist(interactions_nan(ind_triu)./[influence_products(ind_triu)+10^-3], 200)
a = axis; axis([a(1) a(2), 0, sqrt(a(4))]);
xlabel('Synergy normalised by influence product'); ylabel('Count number (close-up)');


% ----
% interactions versus cross-CJ-values

figure(113); clf; set(gca,'FontSize',18)
plot(CV_mean(ind_triu),interactions_nan(ind_triu),'.','Color',[.7 .8 1],'MarkerSize',15);  hold on 
%axis_almost_tight; 
axis tight
a = axis;
plot([a(1) a(2)],[0 0],'--k'); hold on; plot([0 0],[a(3) a(4)],'--k'); hold on;
xlabel('Mean mutual flux control'); ylabel('Synergy');

figure(1113); clf; set(gca,'FontSize',18)
hist(interactions_nan(ind_triu)./[CV_mean(ind_triu)+10^-3], 100)
a = axis; axis([a(1) a(2), 0, sqrt(a(4))]);
xlabel('Synergy normalised by mean mutual flux control'); ylabel('Count number (close-up)');

% ------------------------------------------------------
% histogram of interaction degrees; correlation interaction degrees and viability

n_reactions = length(network.actions)
mean_synergy_degree = mean(synergy_degree) 

figure(114); clf; set(gca,'FontSize',18);
numbers = 0:15;
bar(numbers, [ histc(synergy_degree_pos,numbers) histc(synergy_degree_neg,numbers)],'stacked'); 
legend(sprintf('Positive, total: %s',total_positive_synergies),sprintf('Negative, total: %s',total_negative_synergies)); 
xlabel('Synergy degree'); ylabel('Count number'); 


figure(115); clf; 
plot(influence,synergy_degree,'o'); axis tight; set(gca, 'Fontsize',18);
ylabel('Synergy degree'); xlabel('First-order scaled control coefficient'); 

% sign ratios; 
figure(116); clf 
numbers = -5:5;
bar(numbers, histc(log2_sign_ratio(isfinite(log2_sign_ratio)),numbers)); 
axis tight; set(gca,'FontSize',18);
xlabel('Synergy sign ratio (log2 scale)'); ylabel('Count number');

% -----------------------------------------------------------------------
% interactions between subsystems


if isfield(network,'subsystems'),

  subsystems     = network.subsystems;
  subsystem_list = unique(subsystems);
  ll             = label_names(subsystem_list,subsystems,'multiple');
  ng             = length(network.actions);
  
  for it = 1:length(subsystem_list),
    for it2 = 1:length(subsystem_list),
      M = interaction(ll{it},ll{it2});
      frac_pos(it,it2) = sum(M(:)>0)/prod(size(M));
      frac_neg(it,it2) = sum(M(:)<0)/prod(size(M));
      n_pos(it,it2) = sum(M(:)>0);
      n_neg(it,it2) = sum(M(:)<0);
    end
  end

  % relative numbers of pos / neg interactions between subsystems

  my_M = [n_pos - n_neg]./sqrt([n_pos + n_neg]);

  figure(10000); order_subsystems = sort_by_clustering(my_M); close(10000);

  figure(202); clf; 
  subplot('Position',[0.2 0.02 0.8 0.7]); set(gca,'Fontsize',6);
  im(my_M(order_subsystems,order_subsystems),[],subsystem_list(order_subsystems),subsystem_list(order_subsystems)); 
  colormap(rb_colors);
  colorbar; 
  axis square; axis equal; axis tight; my_xticklabel

  dum = []; dd = 0;
  for it = 1:length(ll),
    dum = [dum; ll{order_subsystems(it)}];
    dd  = [dd; dd(end) + length(ll{order_subsystems(it)})];
  end
  
  figure(201); clf
  im(asinh(interaction(dum,dum))); colormap(rb_colors); colorbar; hold on; 
  for it = 1:length(ll),
    line([0 ng],[dd(it),dd(it)],'Color','k');
    line([dd(it),dd(it)],[0 ng],'Color','k');
  end
  hold off; axis square; axis equal; axis tight
  
end

% -----------------------------------------------------------------------
% hierarchical clustering of (non-thresholded!) interaction profiles

% SCALED EUCLIDEAN YIELDS DIVIDE BY ZERO ERROR!

% which enzymes do have interactions at all?
ind_fin         = find(sum(abs(interaction_thr)));

n_cluster = min(n_cluster,length(ind_fin));
cluster_indices = kmeans(interaction_thr(ind_fin,ind_fin),n_cluster);

% reassign cluster indices such that clusters are order by size
ss = unique(cluster_indices);
clear dd ddd
for it = 1:length(ss), dd(it,1) = sum(cluster_indices==ss(it)); end
[dum,order] = sort(dd);
my_ss = ss(order);
for it = 1:length(order); ddd(cluster_indices==order(it),1) = it; end
cluster_indices = ddd;

% ---------------------------------
% plot on CoHid network

gp = struct('arrowstyle','none', 'linecolor',[0 0 0],'colorbar',0, 'text_offset',[.01,-.01],'colormap',my_spectrum, 'actprintnames',0,'actstyle','fixed','actvaluesmax',max(cluster_indices)/2); % 'squaresize',0.03,

cluster_ind_mapped = nan*zeros(size(network.actions));
cluster_ind_mapped(ind_fin) = cluster_indices;
figure(301); clf; 
netgraph_concentrations(network_CoHid,[], cluster_ind_mapped'-max(cluster_ind_mapped)/2,0,gp); hold off;
set(gcf,'Color',[1 1 1]);

% ---------------------------------
% plot cluster graph NOTE ASINH SCALE

figure(302); clf;
graph_circle_plot(asinh(interaction_thr(ind_fin,ind_fin)),rb_colors,cluster_indices);
axis equal; axis off; 
set(gcf,'Color',[1 1 1]);


% ------------------------------------------------------------
% analysis of cycle signs
% odd cycle lengths  ->  negative signs overrepresented
% even cycle lengths ->  positive signs overrepresented
% (comparison to randomised interaction matrix)

figure(303); clf; 

if flag_analyse_synergy_cycles,

  interaction_neg = interaction_thr; 
  interaction_neg(interaction_neg>0) = 0;

  n_rand = 100;

  [thr_n_cycles,thr_frac_pos,thr_mean_n_cycles,thr_std_n_cycles,thr_pvalue_n_cycles,thr_pvalue_frac_pos] = interaction_count_cycles(interaction_thr,n_rand);
  [neg_n_cycles,neg_frac_pos,neg_mean_n_cycles,neg_std_n_cycles,neg_pvalue_n_cycles] = interaction_count_cycles(interaction_neg,n_rand);
  
  figure(303); clf; 
  subplot(3,1,1); set(gca,'Fontsize',12);
  plot(3:10,thr_n_cycles(3:10),'b*'); hold on;
  errorbar(3:10,thr_mean_n_cycles(3:10),thr_std_n_cycles(3:10),'r*'); hold off;
  set(gca,'YScale','Log'); 
  axis([2.5 10.5 1 10^5]); 
  legend('Cycles from significant synergies','Randomized');
  xlabel('Cycle length'); ylabel('# Cycles'); 
  
  subplot(3,1,2); set(gca,'Fontsize',12);
  plot(3:10,thr_pvalue_frac_pos(3:10),'b*-'); 
  xlabel('Cycle length'); ylabel('P-value for cycles being positive'); 
  
  subplot(3,1,3); set(gca,'Fontsize',12);
  plot(3:10,neg_n_cycles(3:10),'b*'); hold on;
  errorbar(3:10,neg_mean_n_cycles(3:10),neg_std_n_cycles(3:10),'r*'); hold off;
  set(gca,'YScale','Log'); 
  axis([2.5 10.5 1 10^3]); 
  legend('Cycles from significant synergies <0','Randomized');
  xlabel('Cycle length'); ylabel('# Cycles'); 
  
end

% ------------------------------------------------------
% print results

if length(psfile_dir),
  
  display(sprintf('Saving graphics to directory %s',psfile_dir))
  cd(psfile_dir);
  print([es_filenames.psfile_base '_interaction_first_RC.eps'    ], '-f101', '-depsc');
  print([es_filenames.psfile_base '_interaction_second_RC.eps'   ], '-f102', '-depsc');
  print([es_filenames.psfile_base '_interaction_synergy_degree.eps'   ], '-f103', '-depsc');
  print([es_filenames.psfile_base '_interaction_interactions_hist.eps'], '-f111', '-depsc');
  print([es_filenames.psfile_base '_interaction_infl_prod_vs_interaction.eps'], '-f112', '-depsc');
  print([es_filenames.psfile_base '_interaction_CJ_vs_interaction.eps'       ], '-f113', '-depsc');
  print([es_filenames.psfile_base '_interaction_influence_vs_interaction.eps'], '-f1111', '-depsc');
  print([es_filenames.psfile_base '_interaction_infl_prod_vs_interaction_hist.eps'], '-f1112', '-depsc');
  print([es_filenames.psfile_base '_interaction_CJ_vs_interaction_hist.eps'       ], '-f1113', '-depsc');
  print([es_filenames.psfile_base '_interaction_synergy_degree_distribution.eps'], '-f114', '-depsc');
  print([es_filenames.psfile_base '_interaction_influence_vs_synergy_degree.eps'], '-f115', '-depsc');
  print([es_filenames.psfile_base '_interaction_synergy_degree_pos_neg_ratio.eps'], '-f116', '-depsc');

  if isfield(network,'subsystems'),
    print([es_filenames.psfile_base '_interaction_subsystems1.eps'], '-f201', '-depsc');
    print([es_filenames.psfile_base '_interaction_subsystems2.eps'], '-f202', '-depsc');
  end
  
  print([es_filenames.psfile_base '_interaction_clusters.eps'      ], '-f301', '-depsc');
  print([es_filenames.psfile_base '_interaction_cluster_circle.eps'], '-f302', '-depsc');
  if ~flag_analyse_synergy_cycles,
    figure(303)
  end
  print([es_filenames.psfile_base '_interaction_cycle_numbers.eps' ], '-f303', '-depsc');
  
end

figure(101); title('Scaled 1st order resp. coefficients');
figure(102); title('Scaled 2nd order resp. coefficients: fraction(positive)');
figure(103); title('Synergy degrees');
figure(111);  title('Synergies histogram');
figure(1111); title('Influence vs synergy');      
figure(112);  title('Synergy vs  influence product');             
figure(1112); title('Synergy scaled by influence product');      
figure(113);  title('Synergy vs mean CJ'); 
figure(1113); title('Histogram: Synergy / mean CJ ');        
figure(114); title('Distribution of Synergy degrees'); axis tight
figure(115); title('Influence vs Synergy degrees');
figure(116); title('Histogram of synergy log2 sign ratio; result depends strongly on thresholding; here 75 quantile'); 
figure(201); title('Modular Synergies (asinh scale)');
figure(202); title('Fractions of synergy signs (pos-neg)/(pos+neg), Subsystems sorted by clustering');
figure(301); title('Mapped clusters');
figure(302); title('Thresholded synergies (asinh scale) (Reactions without synergies above threshold are omitted)');
figure(303); subplot(3,1,1); title('Cycle numbers (thresholded interaction)');
figure(303); subplot(3,1,2); title('Cycle numbers (aggravating interaction)');
figure(303); subplot(3,1,3); title('P-value (number of positive cycles)');


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other clustering algorithms


% ---------------------------------------------
% compute ordering by eigenvalue decomposition and k-means

neig   = 50;
nclust = 10;

[V,D]                = eig(interaction);
d                    = diag(D);          
[dum,order]          = sort(-abs(d));
X                    = diag(sqrt(abs(d(order)))) * V(:,order(1:neig));
sd                   = sign(d(order(1:neig)));
[cluster_indices, c] = kmeans(X, nclust);
cluster_interactions = c*diag(sd)*c';
[cluster_indices,cluster_interactions] = sort_clusters(network,cluster_indices,cluster_interactions);


% -----------------------------------------------------------------------
% prism algorithm

nclust = 15;

[Z,Conflict,Direct] = prism(interaction_thr,.5,2,1);
cluster_indices     = cluster(Z,'MaxClust',nclust);
cluster_indices     = sort_clusters(network,cluster_indices);


%% OTHER GRAPHICS
% Cumulative distribution of interaction values

% figure(101); clf
% indices = find(triu(ones(size(interaction,1))));
% plot(sort(interaction(indices)),[1:length(indices)]/length(indices));
% %h = bar([-0.05:0.005:0.05],log10(n)); axis tight; 
% %set(h,'FaceColor',[.9 .9 .9]);
% xlabel('Scaled second response coefficients');
% ylabel('Cumulative distribution');

% add titles (after printing)
