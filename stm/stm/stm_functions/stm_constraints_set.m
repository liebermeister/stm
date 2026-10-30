function c = stm_constraints_set(es_constraints,network,varargin)

% STM_CONSTRAINTS_SET Convenience function to define es_constraints
%
% c = stm_constraints_set(es_constraints,network,type1,names1,values1,type2,names2,values2,...)
%
% type:   'v_min', 'v_max', 'v_sign', 'v_fix' (for reactions) or 'ext_signs' (for metabolites)
% names:  list of reaction or metabolite names
% values: vector of values

c = es_constraints;

for it = 1:floor(length(varargin)/3)
  type   = varargin{it*3-2};
  names  = varargin{it*3-1};
  values = varargin{it*3};

  for it2 = 1:length(names)
    switch type
      case {'v_min','v_max','v_sign','v_fix'}, ind = label_names(names(it2),network.actions);
      case {'ext_signs','ext_sign'},           ind = label_names(names(it2),network.metabolites);
      otherwise, error(sprintf('Unknown constraint type %s', type));
    end
    switch type
      case 'v_min',       c.v_min(ind) = values(it2);
      case 'v_max',       c.v_max(ind) = values(it2);
      case 'v_sign',     c.v_sign(ind) = values(it2);
      case 'v_fix',       c.v_fix(ind) = values(it2);
      case {'ext_signs','ext_sign'}, c.ext_signs(ind) = values(it2);
    end
  end

end
