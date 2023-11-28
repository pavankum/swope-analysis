# swope-analysis
Scripts to run Bill Swope's conformer energy analysis
```
NC Number of QM conformers 
NM Number of QM matches with any of the MM structures
NO Number of QM orphans: QM structures that don’t match with any MM struc. within a threshold
EO 10*(num orph)/(num conf): maximum # of confs is 10, might be an equalization metric for # of orphans
PC Num cases match .ne. p-c: Cases where the QM has a close match to a different MM struc (not a parent-child)
LQ 1 if lowest QM is orphan: If the low energy QM does not have any MM matches
LF 1 of lowest FF is orphan: If the MM minima does not correspond to any QM
BF Num FF strucs below ref: Number of FF structures below the reference FF struc that matches to the QM minima
MO FF strucs out of order: number of FF structures out of order rel to QM order
BE Num small dE large ddE: number of cases where ddE is large but dE is small
QS Num 2Q:1F matches: Cases where more than one QM matches to a single FF struc
FS Num 1Q:2F matches: Cases where more than one FF matches to a single QM struc
SC 0-99 score: (# of QM matches - # of flagged cases with BE)/total_qm_conf, higher the better
NL 0-99 score lowE passes: how many have |ddE|<2 among low energy confs (|dE|<4), higher the better 
```

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
