function c = es_constraints_set(es_constraints,network,varargin)

% ES_CONSTRAINTS_SET Convenience function to define es_constraints
%
% c = es_constraints_set(es_constraints,network,varargin)

c = es_constraints;

for it = 1:floor(length(varargin)/3),
  type   = varargin{it*3-2};
  names  = varargin{it*3-1};
  values = varargin{it*3};
  
  for it2 = 1:length(names),
    switch type,
      case {'v_min','v_max','v_sign','v_fix'}, ind = label_names(names(it2),network.actions);
      case 'ext_sign',                         ind = label_names(names(it2),network.metabolites);
    end
    switch type,
      case 'v_min',       c.v_min(ind) = values(it2);
      case 'v_max',       c.v_max(ind) = values(it2);
      case 'v_sign',     c.v_sign(ind) = values(it2);
      case 'v_fix',        c.vfix(ind) = values(it2);
      case 'ext_sign', c.ext_sign(ind) = values(it2);
    end
  end
  
end
