E. coli CCM model: model and data files
=======================================

Modified version of E. coli model (Noor et al 2016)
* Earlier name: e_coli_noor_2016_plus_overflow_oxa_internal.tsv
* Previous version: e_coli_noor_2016_plus_overflow.tsv

Data file e_coli_noor_2016_kinetic_data_modified.tsv
----------------------------------------------------

Derived from the data file e_coli_noor_2016_data.tsv (see github:enzyme-cost-minimization/enzyme-cost-minimization/resource-data/model-files/e_coli_noor_2016)

Same changes that were also made in cellular_enzyme_demands/github/cellular-enzyme-demands/resource-data/models/e_coli_noor_2016_modified/e_coli_noor_2016_ECM_modified_Model.tsv

* Changed PTS / Glucose KM value from 0.0085306826 mM to 0.02 mM (from Millard kinetic E coli model)
  To keep haldane relationship intact: KM value for D_Glucose_6_phosphate 0.20090355 -> 0.47100000
* Changed kcat values for ACN genes:
  substrate catalytic rate constant      ACN_R01325 3.4316827 -> 19.6000
  substrate catalytic rate constant      ACN_R01900 3.3186137 -> 71.0000
  product catalytic rate constant        ACN_R01325 5.2692956 -> 30.0955
  product catalytic rate constant        ACN_R01900 12.502319 -> 267.4806

Data:

* Fluxes are not stationary
* Concentrations are incomplete and have not been tested for thermodynamically feasibility

Notes
-----

* The model is based on the original model (from Noor et al 2016) in
  github:enzyme-cost-minimization/enzyme-cost-minimization/resource-data/model-files/e_coli_noor_2016

* This model was modified, yielding the model version
  github:cellular-enzyme-demands/resource-data/models/e_coli_noor_2016_modified/

* Modifications: 
  o the fbpase reaction was removed
  o Version oxa_internal = 0
    o AcCoA und CoA are set external!!!
    o The pathway is therefore interrupted!
  o Version oxa_internal = 0
    o an overflow reaction (AcCoA -> CoA) was added manually

For usage in structural kinetic modelling, AcetylCoa, CoA, and oxaloacetate have been set internal. A new reaction AcetylCoa <=> CoA
("overflow") and an oxaloacetate-consuming reaction have been introduced to allow for a stationary flux distribution resembling the measured fluxes. 

Afterwards, here, external metabolites Ext1 and Ext2 were added to the new reactions to give flexibility for reaction thermodynamics

Planned future modification: ubiquinone / ubiquinol internal, add oxygen, add a (costly) "oxphos" reaction between them

Note that the parameters stem from the original file e_coli_noor_2016; there is also a variant with updated parameters ("e_coli_noor_2016_ECM_modified_Model.tsv") in the cellular-enzyme-demand project

