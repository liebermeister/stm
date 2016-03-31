function d = es_BASEDIR()

% ES_BASEDIR - Determine directory in which ES toolbox resides
% 
% d = es_BASEDIR()

d = [fileparts(which(mfilename)) '/'];
