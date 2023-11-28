'''
Pavan Behara, Nov 23

Utility script to extract the molecule by industry name tag from the large file
'''
from openeye import oechem
import click

@click.command()
@click.option(
    "--compound_name",
    "compound_name",
    type=click.STRING,
    required=True,
    help="Compound name in industry benchmark set, for example XTP-01128",
)

def main(compound_name):
    ifs = oechem.oemolistream('qm-all-annotated.sdf')
    ofs = oechem.oemolostream(compound_name+'.sdf')
    for mol in ifs.GetOEGraphMols():
         c_name = oechem.OEGetSDData(mol, "Compound Name")
         if compound_name in c_name:
             print(c_name)
             oechem.OEWriteMolecule(ofs, mol)

if __name__=="__main__":
    main()
