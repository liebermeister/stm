function es_dependencies()

% ES_DEPENDENCIES - Check whether the packages used by the STM toolbox are installed
%
% es_dependencies()
%
% (see README.md for details)

if ~exist('mnt_version','file')
  error('Please install the Metabolic Network Toolbox (https://github.com/liebermeister/metabolic-network-toolbox)');
end

% check the packages used by the Metabolic Network Toolbox
% (matlab-utils, SBtab, libSBML, Tensor Toolbox, efmtool)
mnt_dependencies

if ~exist('parameter_balancing','file')
  warning('Please install Parameter balancing for Matlab (https://github.com/liebermeister/parameter-balancing-matlab) - It is used in the workflow demo and for building the E. coli model.');
end
