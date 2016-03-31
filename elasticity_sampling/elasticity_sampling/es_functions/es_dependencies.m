function es_dependencies()

if ~exist('mnt_version','file'),
  error('Please install the Metabolic Network Toolbox (https://github.com/wolframliebermeister/mnt)');
end

mnt_dependencies

if ~exist('tensor','file'),
  error('Please install the Tensor toolbox (see http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html)');
end
