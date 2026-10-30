function es_dependencies(only_required)

% ES_DEPENDENCIES - Check whether the packages used by the STM toolbox are installed
%
% es_dependencies(only_required)
%
% (see README.md for details)
%
% only_required: if true, skip the checks for optional packages (default: false)

if nargin < 1
  only_required = false;
end

if ~exist('mnt_version','file')
  error('Please install the Metabolic Network Toolbox (https://github.com/liebermeister/metabolic-network-toolbox)');
end

if exist('package_version_check','file')
  package_version_check('The Metabolic Network Toolbox', 'mnt_version', '1.3.0');
end

% check the packages used by the Metabolic Network Toolbox
% (matlab-utils, SBtab, libSBML, Tensor Toolbox, efmtool)
mnt_dependencies(only_required)

if ~only_required && ~exist('parameter_balancing','file')
  warning('Please install Parameter balancing for Matlab (https://github.com/liebermeister/parameter-balancing-matlab) - It is used in the workflow demo and for building the E. coli model.');
end

if exist('package_version_check','file') && exist('parameter_balancing_version','file')
  package_version_check('Parameter balancing for Matlab', 'parameter_balancing_version', '1.0.0');
end
