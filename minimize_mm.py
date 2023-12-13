from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange.drivers import get_gromacs_energies, get_openmm_energies
import multiprocessing
import tqdm
import numpy as np
import sys
import openmm
from openff.units.openmm import from_openmm
import os
# suppress stereochemistry warnings
import logging
from rdkit import Chem
import itertools


logging.getLogger("openff").setLevel(logging.ERROR)

def minimize_energy(mol,ff):
    # print(mol)
    # mol.assign_partial_charges('am1bcc',toolkit_registry=AmberToolsToolkitWrapper())
    # mol.assign_partial_charges('am1bccelf10',use_conformers=mol.conformers,toolkit_registry=OpenEyeToolkitWrapper)
    topo = mol.to_topology()
    # print('hi')
    # print(topo.box_vectors)
    interchange = Interchange.from_smirnoff(force_field=ff, topology=topo,allow_nonintegral_charges=True)#,charge_from_molecules=[mol]) #,charge_from_molecules=[mol])
    # interchange.box = unit.Quantity([4,4,4], unit.nanometer)
    integrator = openmm.VerletIntegrator(1 * openmm.unit.femtoseconds)
    simulation = interchange.to_openmm_simulation(integrator)

    # We'll store energies in two lists
    initial_energies = []
    minimized_energies = []

    # And minimized conformers in a second molecule
    minimized_molecule = Molecule.from_topology(topo)
    minimized_molecule.conformers.clear()

    conformer = mol.conformers[0]
    # Tell the OpenMM Simulation the positions of this conformer

    simulation.context.setPositions(conformer.to_openmm())

    # Keep a record of the initial energy
    initial_energies.append(
        simulation.context.getState(getEnergy=True).getPotentialEnergy()
    )

    # Perform the minimization
    simulation.minimizeEnergy()#tolerance=5e-9)

    # Record minimized energy and positions
    min_state = simulation.context.getState(getEnergy=True, getPositions=True)

    minimized_energies.append(min_state.getPotentialEnergy())
    minimized_molecule.add_conformer(from_openmm(min_state.getPositions()))
    return initial_energies,minimized_energies,minimized_molecule

def opt_molecules_parallel(inputs):
    molfile,qmdir,mmdir,ff_file = inputs
    if os.path.isfile('{}/{}'.format(mmdir,molfile)):
        pass
    else:
        sage = ForceField(ff_file)
        mol = Molecule('{}/{}'.format(qmdir,molfile),allow_undefined_stereo=True)
        # mol=Molecule('QM_files/'+molfile,allow_undefined_stereo=True)

        in_e,final_e, final_mol = minimize_energy(mol,sage)
        # print(final_e[0])
        # print(final_e[0].value_in_unit(final_e[0].unit))
        # final_mol.to_file('MM-geometries/{}'.format(molfile),file_format = 'sdf')

        rdkit_final = final_mol.to_rdkit()
        writer = Chem.SDWriter('{}/{}'.format(mmdir,molfile))
        rdkit_final.SetProp('final_energy',str(final_e[0].value_in_unit(openmm.unit.kilocalories_per_mole)))#.value)

        writer.write(rdkit_final)
        # opt_mol.append(final_mol)



if __name__ == '__main__':

    qm_dir='b3lyp-d3bj_dzvp'
    mm_dir = 'openff-2.1.0'
    ff_file_input = 'openff-2.1.0.offxml'
    try: mm_dir=sys.argv[1]
    except IndexError:
        pass
    try: ff_file_input = sys.argv[2]
    except IndexError:
        pass
    try: qm_dir=sys.argv[3] # last since it's not likely to be needed
    except IndexError:
        pass

    if os.path.exists(mm_dir): pass
    else: os.mkdir(mm_dir)

    molfiles = sorted(os.listdir(qm_dir))#[1:2]#[0:1]#.sort()
    # print(molfiles)
    # print(qm_dir)

    with multiprocessing.Pool(8) as pool:
        for x in tqdm.tqdm(pool.imap(opt_molecules_parallel, zip(molfiles,itertools.repeat(qm_dir), itertools.repeat(mm_dir),itertools.repeat(ff_file_input)), chunksize=32),desc = 'Optimizing molecules',total=len(molfiles)):
            pass
