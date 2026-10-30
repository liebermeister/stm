function d = stm_BASEDIR()

% STM_BASEDIR - Determine directory in which ES toolbox resides
%
% d = stm_BASEDIR()

d = [fileparts(which(mfilename)) '/../'];
