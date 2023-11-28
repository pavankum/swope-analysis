'''
Pavan Behara, Nov 23

Script to create the inputs for Swope conformer energies analysis
'''

import json
from openeye import oechem
import pathlib
import click

@click.command()
@click.option(
    "--qm_file_name",
    "qm_file_name",
    type=click.STRING,
    required=True,
    help="File name with the MM optimized conformers, something like qm-all.sdf",
)

@click.option(
    "--mm_file_name",
    "mm_file_name",
    type=click.STRING,
    required=True,
    help="File name with the MM optimized conformers, something like openff-2.0.0-all.sdf",
)
@click.option(
    "--iter_name",
    "iter_name",
    type=click.STRING,
    required=True,
    help="name of the directory will be set like iter_name-data",
)

def main(qm_file_name, mm_file_name, iter_name):
    qm_istream = oechem.oemolistream(qm_file_name)
    mm_istream = oechem.oemolistream(mm_file_name)
    iter = iter_name
    pathlib.Path(f"./{iter}-data/b3lyp-d3bj_dzvp/").mkdir(parents=True, exist_ok=True)
    pathlib.Path(f"./{iter}-data/MM-geometries/").mkdir(parents=True, exist_ok=True)
    compound_list = open(f"./{iter}-data/compound.list", "w")
    
    with open('name_map_to_qca_id.json') as file:
        name_map = json.load(file)
    compound_names = []
    for mol in qm_istream.GetOEGraphMols():
        id = oechem.OEGetSDData(mol, "Record QCArchive")
        name = name_map[id]
        if name[-2:] == '00':
            print(f"Only one conformer for {name}")
            continue
        if name[:-2] not in compound_names:
            compound_names.append(name[:-2])
            compound_list.write(name[:-2]+"\n")
        ofs = oechem.oemolostream(f'{iter}-data/b3lyp-d3bj_dzvp/{name[:-1]}{"%02d"%int(name[-1])}.sdf')
        ofs.SetFormat(oechem.OEFormat_SDF)
        opt_energy = oechem.OEGetSDData(mol, "Energy QCArchive")
        mol.SetTitle(f'{name[:-1]}{"%02d"%int(name[-1])}')
        oechem.OESetSDData(mol, "final_energy", opt_energy)
        oechem.OEWriteMolecule(ofs, mol)
    
    for mol in mm_istream.GetOEGraphMols():
        id = oechem.OEGetSDData(mol, "Record QCArchive")
        name = name_map[id]
        if name[-2:] == '00':
            print(f"Only one conformer for {name}")
            continue
        ofs = oechem.oemolostream(f'{iter}-data/MM-geometries/{name[:-1]}{"%02d"%int(name[-1])}.sdf')
        opt_energy = oechem.OEGetSDData(mol, "Energy FFXML")
        mol.SetTitle(f'{name[:-1]}{"%02d" % int(name[-1])}')
        oechem.OESetSDData(mol, "final_energy", opt_energy)
        ofs.SetFormat(oechem.OEFormat_SDF)
        oechem.OEWriteMolecule(ofs, mol)

if __name__ == "__main__":
    main()

