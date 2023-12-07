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
logging.getLogger("openff").setLevel(logging.ERROR)

def minimize_energy(mol,ff):
    topo = mol.to_topology()
    interchange = Interchange.from_smirnoff(force_field=ff, topology=topo,allow_nonintegral_charges=True) #,charge_from_molecules=[mol])
    # interchange = Interchange.from_smirnoff(force_field=ff, topology=topo)
    interchange.box = unit.Quantity([4, 4, 4], unit.nanometer)
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
    simulation.minimizeEnergy()

    # Record minimized energy and positions
    min_state = simulation.context.getState(getEnergy=True, getPositions=True)

    minimized_energies.append(min_state.getPotentialEnergy())
    minimized_molecule.add_conformer(from_openmm(min_state.getPositions()))
    return initial_energies,minimized_energies,minimized_molecule

def opt_molecules_parallel(molfile):
    if os.path.isfile('MM-geometries/{}'.format(molfile)):
        pass
    else:
        sage = ForceField('openff-2.0.0.offxml')
        mol = Molecule('b3lyp-d3bj_dzvp/'+molfile,allow_undefined_stereo=True)

        in_e,final_e, final_mol = minimize_energy(mol,sage)
        # print(in_e,final_e)
        final_mol.to_file('MM-geometries/{}'.format(molfile),file_format = 'sdf')
        # opt_mol.append(final_mol)

def opt_molecules(mollist):
    print(mollist)
    opt_mol = []
    for molfile in mollist:
        # print(molfile)
        # molname = molfile.split('/')[-1]
        sage = ForceField('openff-2.1.0.offxml')
        mol = Molecule('b3lyp-d3bj_dzvp/'+molfile)

        in_e,final_e, final_mol = minimize_energy(mol,sage)
        # print(in_e,final_e)
        final_mol.to_file('MM-geometries/{}'.format(molfile),file_format = 'sdf')
        opt_mol.append(final_mol)

    return opt_mol

if __name__ == '__main__':
    molfiles = sorted(os.listdir('b3lyp-d3bj_dzvp'))#.sort()
    # molfiles = molfiles[:20]

    with multiprocessing.Pool(8) as pool:
        for x in tqdm.tqdm(pool.imap(opt_molecules_parallel, molfiles),desc = 'Optimizing molecules',total=len(molfiles)):
            pass
