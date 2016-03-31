function r = es_align_models(model_list, field_metabolite_id, field_reaction_id, how)

% [metabolite_ids, reaction_ids, metabolite_indices_list, reaction_indices_list] = es_align_models(model_list,field_metabolite_id,field_reaction_id)
% 
% fields of 'r': [metabolite_ids, reaction_ids, metabolite_indices_list, reaction_indices_list]
% how: 'union' or 'intersection'
% 'intersection' -> items common to all models
% 'union'        -> items appearing in any model

eval(default('how','''intersection''')); %% alternative: 'union'


if 0,
  model_name      = 'eccm'; 
  condition_name  = 'ishii_growth_07';
  filenames       = elasticity_sampling_filenames(model_name,condition_name); 
  cd([filenames.es_dir]); load(filenames.model_data_file);
  model2 = network;
  
  model_name      = 'bscm'; 
  condition_name  = 'glucose';
  filenames       = elasticity_sampling_filenames(model_name,condition_name); 
  cd([filenames.es_dir]); load(filenames.model_data_file);
  model1 = network;
  
  model_name      = 'ycm'; 
  condition_name  = 'diauxic_glucose';
  filenames       = elasticity_sampling_filenames(model_name,condition_name); 
  cd([filenames.es_dir]); load(filenames.model_data_file);
  model3 = network;
  
  model_list = {model1, model2, model3};
  field_metabolite_id = 'metabolite_KEGGID';
  field_reaction_id   = 'reaction_KEGGID';
  
  r = es_align_models(model_list,field_metabolite_id,field_reaction_id);

end


% -----------------------------------------------------------------------------

metabolite_ids   = {};
metabolite_names = {};

for it = 1:length(model_list),
  my_metabolite_ids = getfield(model_list{it},field_metabolite_id);
  my_metabolite_names = model_list{it}.metabolites;
  ind_ok = find(cellfun('length',my_metabolite_ids));
  my_metabolite_ids = my_metabolite_ids(ind_ok);
  my_metabolite_names = my_metabolite_names(ind_ok);
  for itt = 1:length(my_metabolite_ids),
    if ~ismember(my_metabolite_ids{itt},metabolite_ids),
      metabolite_ids    = [metabolite_ids; my_metabolite_ids{itt}];
      metabolite_names  = [metabolite_names; my_metabolite_names{itt}];
    end
  end
end


% -----------------------------------------------------------------------------

reaction_ids   = {};
reaction_names = {};

for it = 1:length(model_list),
  my_reaction_ids = getfield(model_list{it},field_reaction_id);
  my_reaction_names = model_list{it}.actions;
  ind_ok = find(cellfun('length',my_reaction_ids));
  my_reaction_ids   = my_reaction_ids(ind_ok);
  my_reaction_names = my_reaction_names(ind_ok);
  for itt = 1:length(my_reaction_ids),
    if ~ismember(my_reaction_ids{itt},reaction_ids),
      reaction_ids    = [reaction_ids; my_reaction_ids{itt}];
      reaction_names  = [reaction_names; my_reaction_names{itt}];
    end
  end
end

% -----------------------------------------------------------------------------

for it = 1:length(model_list),
  my_metabolite_ids = getfield(model_list{it},field_metabolite_id);
  metabolite_indices_list(:,it) = label_names(metabolite_ids , my_metabolite_ids);
end

for it = 1:length(model_list),
  my_reaction_ids = getfield(model_list{it},field_reaction_id);
  reaction_indices_list(:,it) = label_names(reaction_ids , my_reaction_ids);
end

% -----------------------------------------------------------------------------

switch how, 
  case 'union',
  case 'intersection',
    ind_ok = find(prod(metabolite_indices_list,2));
    metabolite_ids = metabolite_ids(ind_ok);
    metabolite_names = metabolite_names(ind_ok);
    metabolite_indices_list = metabolite_indices_list(ind_ok,:);
    ind_ok = find(prod(reaction_indices_list,2));
    reaction_ids = reaction_ids(ind_ok);
    reaction_names = reaction_names(ind_ok);
    reaction_indices_list = reaction_indices_list(ind_ok,:);
end


% -----------------------------------------------------------------------------

r.metabolite_ids          = metabolite_ids          ;  
r.metabolite_names        = metabolite_names        ;
r.reaction_ids            = reaction_ids            ;
r.reaction_names          = reaction_names          ;
r.metabolite_indices_list = metabolite_indices_list ;
r.reaction_indices_list   = reaction_indices_list   ;