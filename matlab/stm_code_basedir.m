function d = stm_code_basedir()

% d = stm_code_basedir()
% returns the path to this directory (used for file paths)
  
d = [fileparts(which(mfilename))];
