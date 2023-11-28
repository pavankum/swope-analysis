# swope-analysis
Scripts to run Bill Swope's conformer energy analysis

# File manifest
 create_inputs_for_swope_analysis.py - takes in the QM and MM files (both ordered the same, n^{th} molecule in each are the same)

 name_map_to_qca_id.json - name map from QCA id to industry benchmark set compound id

 computermsd.f - computes the QM x MM RMSD matrix

 tfd_matrix_qm_by_mm.py - creates a QM x MM matrix in a format similar to the fortran script which calculates RMSD

 computematch.f - computes the fingerprint based on the QM x MM RMSD/TFD matrix for each molecule

 computematch.out - compiled executable

 computermsd.out - compiled executable

 get_mol_with_name.py - utility script to extract a molecule by name

 MM-all.sdf - MM optimized conformer geometries in SDF format

 qm-all-annotated.sdf - QM reference conformer geometries in SDF format

 ./test-data - directory created with the create_inputs_for_swope_analysis.py script



# compile the fortran files
Compile the fortran files `computermsd.f` and `computematch.f` using gfortran.
```
gfortran computermsd.f -llapack -o computermsd.out
gfortran computematch.f -llapack -o computematch.out
```

# Scripts to create the input files from the QM and MM files
QM references in one file and corresponding MM reference in the other. This script will create directories in the format needed by the above scripts (compute\*)
```
python create_inputs_for_swope_analysis.py --qm_file_name qm-all-annotated.sdf --mm_file_name MM-all.sdf --iter_name test
```

# run the bash script to generate report
Copy the bash script `report.sh` to the data folder and run it
