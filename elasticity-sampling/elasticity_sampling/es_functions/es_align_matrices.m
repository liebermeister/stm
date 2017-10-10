function aligned_matrices = es_align_matrices(matrix_list,r,type_rows,type_columns)

%ES_ALIGN_MATRICES Align matrices corresponding to different metabolic models
%
% aligned_matrices = es_align_matrices(matrix_list,r,type_rows,type_columns)
% 
% fields of 'r': [metabolite_ids, reaction_ids, metabolite_indices_list, reaction_indices_list]
% can be prepared by es_align_models
%
% type_rows, type_columns: either 'metabolite' or 'reaction'

if 0,
  %% example: align the stoichiometric matrices of three models
  compare_eccm_bscm_ycm
end

nm = length(r.metabolite_ids);
nr = length(r.reaction_ids);

switch type_rows,
  case 'metabolite', n_rows = nm; ind_rows = r.metabolite_indices_list;
  case 'reaction',   n_rows = nr; ind_rows = r.reaction_indices_list;
  case 'none',       n_rows = 1;  ind_rows = ones(1,length(matrix_list));
end

switch type_columns,
  case 'metabolite', n_columns = nm; ind_columns = r.metabolite_indices_list; 
  case 'reaction',   n_columns = nr; ind_columns = r.reaction_indices_list; 
  case 'none',       n_columns = 1;  ind_columns = ones(1,length(matrix_list));
end


for it = 1:length(matrix_list),
  aligned_matrices{it} = nan * zeros(n_rows, n_columns);
  my_ind_rows    = ind_rows(:,it);
  my_ind_columns = ind_columns(:,it);
  matrix_list{it}(my_ind_rows(find(my_ind_rows)),my_ind_columns(find(my_ind_columns)));
  aligned_matrices{it}(find(my_ind_rows),find(my_ind_columns)) = matrix_list{it}(my_ind_rows(find(my_ind_rows)),my_ind_columns(find(my_ind_columns)));
end
