function [ext_sign_vector,v_sign_vector] = make_ext_sign_vector_vector(ext_sign,metabolites,ext_ind,v_signs,actions)

%MAKE_EXT_SIGN_VECTOR_VECTOR - Helper function
%
% [ext_sign_vector,v_sign_vector] = make_ext_sign_vector_vector(ext_sign,metabolites,ext_ind,v_signs,actions)

nm              = length(metabolites);
ext_sign_vector = nan*ones(nm,1);

for it=1:length(ext_sign),
  this_ind = find(strcmp(metabolites,ext_sign{it}{1}));
  if find(this_ind==ext_ind),
    ext_sign_vector(this_ind) = ext_sign{it}{2};
  end
end

if nargout >1,

  nr            = length(actions);
  v_sign_vector = nan*ones(nr,1);

  for it=1:length(v_signs),
    this_ind = find(strcmp(actions,v_signs{it}{1}));
    v_sign_vector(this_ind) = v_signs{it}{2};
  end
  
end