## Elastic and Anelastic Gaussian Perturbations at Single Lat, Lon, Depth and Given Sigma

### Instructions

1. Compile SPECFEM3D_GLOBE with ADIOS
2. Run Mesher by turning ADIOS one. The solver_data.bp file will be used from the DATABASES_MPI directory
3. Compile code by changing the SPECFEM directory, ADIOS directory and also the add the OUTPUT_FILES directory of your SPECFEM run in the $SPECFEM_INC of Makefile.
4. make all 
5. Usage Instruction:
   First create the perturbed Gaussian Model saved in *model_pert.bp* of the given directory
   *Usage: bin/xmodel_pert lat lon depth sigma solver_dir output_dir 'array specifying which models to perturb in the order vph vpv vsh vsv rho eta qmu'*


   Create only perturbation for elastic structure, Model file saved as *gll_model_perturbed.bp*
   *Usage: bin/xmodel_add model_new_directory model_old_directory model_perturbation_directory MaximumPerturbation*


   
   Create perturbation in anelastic structure, need to run elastic perturbation before running this, output/new directory of xmodel_add should be the input/old model for xmodel_add_qmu, Final Model is saved as *model_gll_perturbed_qmu.bp*
   *Usage: bin/xmodel_add_qmu model_new_directory model_old_directory model_perturbation_directory MaximumPerturbation*
 
