-------------------------------
Glycolysis in human hepatocytes
-------------------------------


Model construction: 
-------------------

1. The original kinetic model (MATLAB file) was obtained from the author of [1].

2. Starting from the original substance concentrations, a steady state was computed, 

3. A new model was constructed, based on the original model's network structure 
   and steady state, by elasticity sampling (data in struct 'reference_state' in files .._reference_state). 


Model files:
------------

1. Hepatocyte_Glycolysis_ES.xml                  SBML file
2. Hepatocyte_Glycolysis_ES.mat                  MATLAB .mat file (for MNT toolbox)
3. Hepatocyte_Glycolysis_ES_reference_state.mat  MATLAB .mat file (for MNT toolbox)


[1] M. König and H.G. Holzhütter (2012)
    Kinetic modeling of human hepatic glucose metabolism in T2DM predicts 
    higher risk of hypoglycemic events in rigorous insulin therapy. 
    J Biol Chem., 287:36978-36989, 2012.
