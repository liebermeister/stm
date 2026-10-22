function setup_test_path()

% setup_test_path()
%
% Add the toolbox to the Matlab path and check that the required packages are installed
% (called at the start of each test file)

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'matlab')));

if ~exist('mnt_version','file')
  error('Please add the Metabolic Network Toolbox (https://github.com/liebermeister/metabolic-network-toolbox) to the Matlab path');
end

if ~exist('lp_solve','file')
  error('Please add matlab-utils (https://github.com/liebermeister/matlab-utils) to the Matlab path');
end
