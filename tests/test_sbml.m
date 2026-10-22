% Tests for building a model from an SBML file (requires libSBML)

setup_test_path;
sbml_file = fullfile(fileparts(mfilename('fullpath')), '..', 'resource-data', 'models-article', ...
                     'Example_Linear_Pathway', 'Example_Linear_Pathway_model.xml');

%% es_sbml_to_model
if ~exist('TranslateSBML', 'file')
  warning('libSBML is not installed; skipping the test');
  return
end
[network, result] = es_sbml_to_model(sbml_file);
assert(isfield(network, 'kinetics'));
N = network.N(network.external == 0, :);
assert(max(abs(N * result.v)) < 1e-8);
