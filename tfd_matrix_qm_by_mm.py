import json
import os
from openeye import oechem
import pathlib
import click
from rdkit.Chem import TorsionFingerprints
from openff.toolkit.topology import Molecule
import numpy as np

@click.command()
@click.option(
    "--dir_name",
    "dir_name",
    type=click.STRING,
    required=True,
    help="name of the directory to look into for inputs and create a TFD directory",
)

# The input directories are hard-coded to match Bill Swope's analysis scripts
# Passing the directory name that already has the QM and MM structures within dir_name/b3lyp-d3bj_dzvp, dir_name/openff-2.0.0, 
# and acts on a list of compounds in dir_name/compound.list
# This creates a directory TFD-OpenFF with compund_name.tfd files that are similar in structure to rmsd files
# The matrix is QM x MM TFDs 

def main(dir_name):
    pathlib.Path(f"{dir_name}/TFD-OpenFF/").mkdir(parents=True, exist_ok=True)
    with open(dir_name+'/compound.list', 'r') as f:
        compounds = f.read().splitlines() 
    files = [filename for filename in os.listdir(dir_name+'/b3lyp-d3bj_dzvp/') if os.path.isfile(os.path.join(dir_name+'/b3lyp-d3bj_dzvp/',filename))]
    for compound_name in compounds:
        prefixed = [filename for filename in files if filename.startswith(compound_name)]
        prefixed = sorted(prefixed)
        tfds = np.empty((len(prefixed), len(prefixed)))
        tfds.fill(np.nan)
        # read molecule files
        openff_qm_mols = [Molecule.from_file(dir_name+'/b3lyp-d3bj_dzvp/'+file, file_format='sdf', allow_undefined_stereo=True) for file in prefixed]
        openff_mm_mols = [Molecule.from_file(dir_name+'/openff-2.0.0/'+file, file_format='sdf', allow_undefined_stereo=True) for file in prefixed]
        # convert to rdmols
        rd_qm_mols = [openff_qm_mol.to_rdkit() for openff_qm_mol in openff_qm_mols]
        rd_mm_mols = [openff_mm_mol.to_rdkit() for openff_mm_mol in openff_mm_mols]
        # energies
        qm_energies = [openff_qm_mol.properties['Energy QCArchive'] for openff_qm_mol in openff_qm_mols]
        mm_energies = [openff_mm_mol.properties['Energy FFXML'] for openff_mm_mol in openff_mm_mols]  
        skip_record = False
        for i, rd_qm_mol in enumerate(rd_qm_mols):
            for j, rd_mm_mol in enumerate(rd_mm_mols):
                try:
                    tfd = TorsionFingerprints.GetTFDBetweenMolecules(
                        rd_qm_mol, rd_mm_mol
                    )

                except IndexError:

                    print(
                        f"error calculating TFD for id={record_id} - possibly no "
                        f"non-terminal rotatable bonds found."
                    )
                    tfd = np.nan
                except ValueError:
                    print("Mismatch in conformer geometeries, rdkit thinks they are different molecules\n")
                    print("Failed for " + compound_name)
                    tfd = np.nan
                    skip_record = True
                    break


                tfds[i, j] = tfd
        if not skip_record:
            with open(dir_name+'/TFD-OpenFF/'+compound_name+'.tfd', 'w') as f:
                f.write(str(len(prefixed))+'\n')
                for i, item in enumerate(prefixed):
                    f.write(dir_name+'/b3lyp-d3bj_dzvp/'+item+'\n')
                    f.write('\t'+qm_energies[i]+'\n')
                f.write(str(len(prefixed))+'\n')
                for i, item in enumerate(prefixed):
                    f.write(dir_name+'/openff-2.0.0/'+item+'\n')
                    f.write('\t'+mm_energies[i]+'\n')
                np.savetxt(f, tfds, fmt='%-2.4f')


if __name__ == "__main__":
    main()

