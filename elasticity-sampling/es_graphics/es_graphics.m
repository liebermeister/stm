function M = es_graphics(network, network_CoHid, graphics_type, data, p)

% M = es_graphics(network, network_CoHid, graphics_type, data, p)
%
% graphics_type: names, fluxes, flux_movie, ...

eval(default('data','struct','p','struct'));

p_default = struct('print',0,'show_regulation',0, 'metprintnames', 0);
p         = join_struct(p_default, p);

if p.print, 
  %network_CoHid.graphics_par.FontSize = 6;
  %network.graphics_par.FontSize       = 6;
end

% fix (to be removed later)
network_CoHid.graphics_par = join_struct(struct('omitreactions',[],'omitmetabolites',[]),network_CoHid.graphics_par);

ii = setdiff(1:length(network_CoHid.actions),label_names(network_CoHid.graphics_par.omitreactions,network_CoHid.actions));
ind_reactions_show   = network_CoHid.graphics_par.reaction_mapping(ii);

ii = setdiff(1:length(network_CoHid.metabolites),label_names(network_CoHid.graphics_par.omitmetabolites,network_CoHid.metabolites));
ind_metabolites_show = network_CoHid.graphics_par.metabolite_mapping(ii);

gp = join_struct(struct('metstyle','fixed','actprintnames',0,'arrowstyle','fluxes','colormap',rb_colors,'colorbar',1,'split_names',0,'arrowcolor',[0.7, 0.7, 0.7],'show_regulation',p.show_regulation, 'FontSize', 8, 'colorbar_fontsize',20),p);

if isfield(data,'cb_numbers'), 
  gp.colorbar = 1; 
  gp.colorbar_numbers = data.cb_numbers; 
end 

M = [];

switch graphics_type,
  
  case 'names',
    gp.metprintnames = 1; 
    gp.squaresize    = 0.5 * network_CoHid.graphics_par.squaresize;
    gp.actprintnames = 0; 
    gp.arrowvalues   = [];% ones(size(network.actions)); 
    gp.arrowstyle    = 'directions';
    gp.colorbar      = 0;
    netgraph_concentrations(network_CoHid,ones(size(network.metabolites)),[],1,gp);

  case 'metabolites',
    gp.metprintnames    = p.metprintnames; 
    gp.metvaluesmax     = log10(100);
    %gp.colorbar         = 0;
    gp.colorbar_numbers = log10([0.01,0.1,1,10,100]);
    netgraph_concentrations(network_CoHid,log10(data.c),[],1,gp);

  case 'metabolite_ratios',
    gp.metprintnames    = p.metprintnames; 
    gp.metvaluesmax     = log10(100);
    gp.colorbar_numbers = [0.01,0.1,1,10,100];
    netgraph_concentrations(network_CoHid,log10(data.c),[],1,gp);
  
  case 'enzymes',
    gp.metprintnames    = p.metprintnames; 
    gp.showsign         = 0;
    gp.actvaluesmin     = log10(0.000001);
    gp.actvaluesmax     = log10(1);
    gp.colorbar_numbers = [0.000001 0.001 1];
    gp.arrowstyle       = 'none';
    gp.actstyle         = 'fixed';
    netgraph_concentrations(network_CoHid,[],log10(data.u),1, gp);
  
  case 'enzyme_ratios',
    gp.metprintnames    = p.metprintnames; 
    gp.showsign         = 0;
    gp.actvaluesmax     = log10(100);
    gp.colorbar_numbers = [0.01,0.1,1,10,100];
    gp.arrowstyle = 'none';
    gp.actstyle = 'fixed';
    gp.colorbar         = 1;
    netgraph_concentrations(network_CoHid,[],log10(data.u),1, gp);
  
  case 'fluxes',
    gp.metprintnames   = p.metprintnames; 
    gp.arrowvalues     = data.v / max(abs(data.v));%max(abs(data.v(ind_reactions_show)));
    gp.arrowvalues(isnan(gp.arrowvalues))= 0;
    gp.single_arrow    = 1;
    gp.arrow_stoichiometries = 0;
    gp.actstyle        = 'none';
    gp.colorbar        = 1;
    gp.arrowcolor      = [1 0 0];
    %% show SIGNS OF production rates
    if isfield(data,'production_rates'),
      data.production_rates = sign(data.production_rates) .* double(abs(data.production_rates)>10^-8);
      max_prod_rate   = max(abs(data.production_rates(ind_metabolites_show)));
      gp.metvaluesmax = max(max_prod_rate, 10^-10);
      netgraph_concentrations(network_CoHid,data.production_rates,[],1,gp);%data.v
    else,
      gp.colorbar = 0;
      netgraph_concentrations(network_CoHid,[],[],1,gp);%abs(data.v); %% -zero_to_nan(network.external)
    end

  case 'elasticities',
    gp.metprintnames  = p.metprintnames; 
    gp.squaresize  = 0.8 * network_CoHid.graphics_par.squaresize;
    gp.arrowvalues = sign(data.v);
    gp.arrowvaluesmax = 1; % max(abs(data.v(ind_reactions_show)));
    gp.arrowsize  = 0.5 * network_CoHid.graphics_par.arrowsize;
    gp.arrowstyle  = 'fluxes';
    gp.edgevalues  = data.elasticities';
    gp.edgestyle   = 'fixed';
    gp.colorbar_numbers = max(abs(full(data.elasticities(:)))) * [-1,0,1];%[min(data.elasticities(:)) max(data.elasticities(:))];
    gp.show_regulation  = 1;
    gp.show_regulationvalues = 0;
    gp.regulationvalues = data.elasticities;
    gp.regulationstyle = 'fixed';
    netgraph_concentrations(network_CoHid,[],[],1,gp);    

  case 'flux_movie',
    gp.metprintnames    = 1; 
    M = netgraph_flux_movie(network_CoHid,zero_to_nan(network.N*data.v,10^-5),data.v,1,data.n_frames);
  
  case 'thermodynamic_forces',
    gp.metprintnames    = p.metprintnames; 
    theta_directed = data.theta; 
    theta_directed(data.v<0) = - theta_directed(data.v<0);
    theta_directed(data.v==0) = nan;
    gp.colorbar_numbers = [min(theta_directed) max(theta_directed)];
    gp.arrowstyle       = 'none';
    gp.actstyle         = 'box';
    %gp.colorbar        = 0;
    netgraph_concentrations(network_CoHid,[],theta_directed,1,gp);

  case 'dissipation',
    gp.metprintnames    = p.metprintnames; 
    sigma = RT * data.theta .* data.v; 
    gp.colorbar_numbers = [min(sigma) max(sigma)];
    gp.arrowstyle       = 'none';
    %gp.colorbar         = 0;
    netgraph_concentrations(network_CoHid,[],sigma,1,gp);

  case 'chemical_potentials',
    gp.metprintnames    = p.metprintnames; 
    gp.metvaluesmin = min(data.mu(ind_metabolites_show)); 
    gp.metvaluesmax = max(data.mu(ind_metabolites_show));
    gp.colorbar_numbers = [min(data.mu) max(data.mu)];
    %gp.colorbar         = 0;
    gp.showsign         = 0;
    netgraph_concentrations(network_CoHid,data.mu,[],1,gp);    
    
  case 'response_coefficients',
    gp.metprintnames    = p.metprintnames; 
    gp = struct('actvalues', data.R, 'actstyle','fixed','arrowvalues',data.v,'arrowcolor',[0.7,0.7,0.7], 'colorbar', 1, 'colormap',rb_colors);
    gp.arrowsize      = network_CoHid.graphics_par.arrowsize;
    gp.arrowvalues    = data.v / max(abs(data.v(ind_reactions_show)));
    %% gp.arrowstyle  = 'none';
    %% gp.arrowsize  = 0.5 * network_CoHid.graphics_par.arrowsize;
    %% gp.arrowvalues = data.v;
    %% gp.arrowvaluesmax   = max(abs(data.v(ind_reactions_show)));
    %%    dum = sort(data.R); gp.actvaluesmax = dum(end);
    gp.actvaluesmax     = max(abs(data.R(ind_reactions_show)));
    gp.colorbar_numbers = [-1 0 1];
    netgraph_concentrations(network_CoHid,[],data.R,1,gp);
  
  case 'response_coefficients_Sext',
    gp.metprintnames  = p.metprintnames; 
    gp                = struct('metstyle','fixed','arrowcolor',[0.7,0.7,0.7], 'colorbar',1, 'colormap',rb_colors);
    gp.arrowsize      = network_CoHid.graphics_par.arrowsize;
    gp.arrowvalues    = data.v / max(abs(data.v(ind_reactions_show)));
    %% gp.arrowstyle  = 'none';
    %% gp.arrowsize  = 0.5 * network_CoHid.graphics_par.arrowsize;
    %% gp.arrowvalues = data.v;
    %% gp.arrowvaluesmax   = max(abs(data.v(ind_reactions_show)));
    gp.metvaluesmax = max(abs(data.R(ind_metabolites_show)));
    netgraph_concentrations(network_CoHid,zero_to_nan(data.R),[],1,gp);
    
  case 'individual_metabolites',
    if isfield(data,'log_c'),
      ind_data = find(sum(isfinite(data.log_c.mean),2));
      n_exp = size(data.log_c.mean,2);
    else,
      ind_data = 1:size(network.N,1);
      n_exp = size(data.model(1).c,2);      
    end
    [ni,nk] = subplot_n(length(ind_data));
    for it = 1:length(ind_data),
      this_ind = ind_data(it);
      subplot(ni,nk,it); 
      hold on;
      if isfield(data,'model'),
        for it = 1: prod(size(data.model)),
          plot(1:n_exp,log(data.model(it).c(this_ind,:))/log(10),'c-');
        end
      end
      if isfield(data,'log_c'),
        errorbar(1:n_exp,data.log_c.mean(this_ind,:)/log(10),data.log_c.std(this_ind,:)/log(10),'ro');
      end
      hold off
      title(network.metabolites{this_ind}); axis tight; %ylabel('concentration (log_{10} mM)'); %
    end
    
  case 'individual_fluxes',
    if isfield(data,'genes'), names = network.genes; unknown = find(strcmp(names,'nan')); names(unknown)= network.actions(unknown);
    else, names = network.actions; end
    ind_data = find(sum(isfinite(data.v.mean),2));
    n_exp = size(data.v.mean,2);
    [ni,nk] = subplot_n(length(ind_data));
    for it = 1:length(ind_data),
      this_ind = ind_data(it);
      subplot(ni,nk,it);         hold on;
      if isfield(data,'model'),
        for it = 1: prod(size(data.model)),
          plot(1:n_exp,data.model(it).v(this_ind,:),'c-');
        end
      end
      errorbar(1:n_exp,data.v.mean(this_ind,:),data.v.std(this_ind,:),'ro'); 
      hold off
      title(names{this_ind}); ylabel('velocity (mM/s)'); axis tight; 
    end
  
  case 'individual_enzymes',
    if isfield(data,'genes'), names = network.genes; unknown = find(strcmp(names,'nan')); names(unknown)= network.actions(unknown);
    else, names = network.actions; end
    ind_data = find(sum(isfinite(data.log_u.mean),2));
    n_exp = size(data.log_u.mean,2);
    [ni,nk] = subplot_n(length(ind_data));
    for it = 1:length(ind_data),
      this_ind = ind_data(it);
      subplot(ni,nk,it); 
        hold on;
        if isfield(data,'model'),
          for it = 1: prod(size(data.model)),
            plot(1:n_exp,log(data.model(it).u(this_ind,:))/log(10),'c-');
          end
        end
        errorbar(1:n_exp,data.log_u.mean(this_ind,:)/log(10),data.log_u.std(this_ind,:)/log(10),'ro'); 
        hold off
        title(names{this_ind}); ylabel('enzyme concentration (log_{10} mM)'); axis tight;
    end
end
