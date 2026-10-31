function d = stm_BASEDIR()

% STM_BASEDIR - Determine directory in which the STM toolbox resides
%
% d = stm_BASEDIR()

d = [fileparts(which(mfilename)) '/../'];
