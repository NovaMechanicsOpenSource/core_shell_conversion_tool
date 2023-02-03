
'''This is a python script to add or remove shell particles from a data/restart LAMMPS file.
It is heavily based on the LAMMPS python API so that the required code can be kept to a minimum.'''

import sys
import os
from argparse import ArgumentParser
import numpy as np
from lammps import lammps
from particle_system import ParticleSystem

def create_cli():
    '''Create an elementary CLI based on the argparse module & perform initial consistency checks'''

    # set the command line interface
    parser = ArgumentParser(description='Removing/adding shell particles to a LAMMPS system')

    parser.add_argument("input_file", help="The full path & name of the input data/restart file")
    parser.add_argument("output_file", help="The full path & name of the output data file")
    parser.add_argument("shell_atoms", nargs='+', default=[], type=int, help="The atom type\
            of the shell particles which will be removed")
    parser.add_argument("-r","--restart_file", action='store_true', help="Indicates if the\
            input file is restart or not")
    parser.add_argument("-e","--extra_keywords", default=[], help="Extra keywords to be\
            used when reading data from LAMMPS input file")
    parser.add_argument('-a',"--atom_style", default="full", help="Specifies the atom_style\
            of the input/output files")
    parser.add_argument("-f","--forcefield", default=None, help="The full path & name of the\
            force field file if shell particles are added")

    args = parser.parse_args()

    # check that the specified input file exists
    if not os.path.exists(args.input_file):
        print("The specified LAMMPS input file does not exist. The script will exit")
        sys.exit()

    # check if the specified output file exists and issue a warning
    if os.path.exists(args.output_file):
        print("The specified LAMMPS output file already exists. The script will overwrite it")

    # check that the specified atom style supports partial charges and bonds
    dummy_lmp = lammps()
    dummy_lmp.command(f"atom_style {args.atom_style}")

    if dummy_lmp.extract_global("q_flag") == 0: # case for partial charge
        print("The specified atom style does not support partial charges. The script will exit")
        sys.exit()

    if dummy_lmp.extract_global("molecule_flag") == 0: # case for bond
        print("The specified atom style does not support bonds. The script will exit")
        sys.exit()

    dummy_lmp.close()

    if args.forcefield is not None:
        if not os.path.exists(args.forcefield):
            print("The specified force field file does not exist. The script will exit")
            sys.exit()

    return args

def polarize(cur_sys, args, ff_atom_type, ff_mass_fraction, ff_partial_charge):
    '''Create the shell particles'''

    n_new_types = len(args.shell_atoms)

    # define a small number
    small_number = np.finfo(float).eps

    # split the masses of the atoms
    offset_for_types = cur_sys.ntypes - n_new_types
    for jtype, itype in enumerate(args.shell_atoms):
        imass = cur_sys.mass_per_type[itype-1]
        fraction = ff_mass_fraction[ff_atom_type == itype][0]
        cur_sys.mass_per_type[itype - 1] = fraction*imass
        cur_sys.mass_per_type[jtype + offset_for_types] = (1 - fraction)*imass

    # create the shell particles for each 'full' atom
    existing_atoms = cur_sys.natoms

    for j, itype in enumerate(args.shell_atoms):
        indx = (cur_sys.atom_type == itype)
        n_atoms = np.count_nonzero(indx)
        n_atom_id = [existing_atoms+iat+1 for iat in range(0, n_atoms)]
        n_atom_type = [j+offset_for_types+1] * n_atoms
        n_atom_coords = [0]*(3*n_atoms)
        n_atom_coords[0::3] = cur_sys.atom_x_coordinates[indx] + \
                4*small_number*(np.random.random_sample((n_atoms,)) - 0.5)
        n_atom_coords[1::3] = cur_sys.atom_y_coordinates[indx] + \
                4*small_number*(np.random.random_sample((n_atoms,)) - 0.5)
        n_atom_coords[2::3] = cur_sys.atom_z_coordinates[indx] + \
                4*small_number*(np.random.random_sample((n_atoms,)) - 0.5)
        existing_atoms += n_atoms

        new_atoms = cur_sys.lmp_pointer.create_atoms(n_atoms, n_atom_id, n_atom_type, n_atom_coords)
        if new_atoms != n_atoms:
            print("There is a problem in creating the new atoms ! The script will exit")
            sys.exit()

    # create the bonds between the core and shell particles
    btype = []
    batom1 = []

    nbondtypes = cur_sys.nbondtypes - n_new_types

    # populate the lists with the bond type & the pair of core-shell particles
    for j, itype in enumerate(args.shell_atoms):
        btype.extend(nbondtypes+j+1  for _ in cur_sys.atom_id[cur_sys.atom_type == itype])
        batom1.extend(cur_sys.atom_id[cur_sys.atom_type == itype])

    batom2 = list(range(cur_sys.natoms+1, existing_atoms+1))

    # send the information to lammps using the create_bonds command
    for ibtype, ibatom1, ibatom2 in zip(btype, batom1, batom2):
        cur_sys.lmp_pointer.command(f"create_bonds single/bond {ibtype} {ibatom1} {ibatom2} \
                special no")

    # fix the partial charges of the atoms
    # set the partial charges of all new atom types to zero
    cur_sys.lmp_pointer.command(f"set type {offset_for_types+1}* charge 0.0")

    # extend the arrays for storing the atomic properties without calling LAMMPS
    cur_sys.atom_partial_charge.resize(existing_atoms)
    cur_sys.atom_type.resize(existing_atoms)
    cur_sys.atom_mol.resize(existing_atoms)
    cur_sys.atom_id.resize(existing_atoms)

    # update the ids
    cur_sys.atom_id[cur_sys.natoms:] = batom2

    # update the atom type
    existing_atoms = cur_sys.natoms
    for j, itype in enumerate(args.shell_atoms):
        cur_sys.atom_type[existing_atoms:existing_atoms+n_atoms] = j+offset_for_types+1
        existing_atoms += n_atoms

    # correct the partial charge
    for indx, itype in enumerate(args.shell_atoms):
        q_incr = ff_partial_charge[ff_atom_type == itype][0]
        cur_sys.lmp_pointer.command(f"set type {offset_for_types + indx + 1} charge {q_incr}")

    for iat in batom1:
        idx = np.where(cur_sys.atom_id == iat)[0][0]
        itype = cur_sys.atom_type[idx]
        q_incr = ff_partial_charge[ff_atom_type == itype][0]
        cur_sys.lmp_pointer.command(f"set atom {iat} charge \
                {cur_sys.atom_partial_charge[idx]-q_incr}")

    # correct the molecule property of each atom
    if cur_sys.atom_mol is not None:
        for iat, jat in zip(batom1, batom2):
            idx = np.where(cur_sys.atom_id == iat)[0][0]
            cur_sys.lmp_pointer.command(f"set atom {jat} mol {cur_sys.atom_mol[idx]}")

    # send the information for the mass per type to lammps
    cur_sys.set_mass_per_type()

    # print some information for the core-shell treatment by the force field
    print(f"The types of the created shell particles are \
            {' '.join([str(x) for x in range(offset_for_types+1, cur_sys.ntypes + 1)])}")
    print(f"The bond types for the core-shell oscillators are \
            {' '.join([str(x) for x in range(nbondtypes+1, nbondtypes+n_new_types + 1)])}")

def core_shell():
    '''Adding/Removing shell particles from a lammps data/restart file'''

    # set the command line interface & perform initial consistency checks
    args = create_cli()

    # Allocate extra space if shell particles will be create & read the force field file
    n_new_types = 0
    if args.forcefield is not None:
        n_new_types = len(args.shell_atoms)

    # Initialization
    lmp = lammps()
    lmp.command(f"atom_style {args.atom_style}") # define the atom_style to use

    # Read the input file
    extra_keyword = ' '.join([str(x) for x in args.extra_keywords])
    if args.restart_file:
        # Restart input file
        lmp.command(f"read_restart {args.input_file} extra/atom/types {n_new_types} \
                extra/bond/types {n_new_types} extra/bond/per/atom 1 nocoeff {extra_keyword} ")
    else:
        # Data input file
        lmp.command(f"read_data {args.input_file} extra/atom/types {n_new_types} \
                extra/bond/types {n_new_types} extra/bond/per/atom 1 nocoeff {extra_keyword} ")

    # set the communication cutoff to a non-zero value
    box_info = lmp.extract_box()
    lmp.command(f"pair_style zero {0.1*np.min(np.asarray(box_info[1])-np.asarray(box_info[0]))} ")
    lmp.command("bond_style zero nocoeff")
    lmp.command("pair_coeff * *")
    lmp.command("bond_coeff *")

    # define the object with all system data
    current_system = ParticleSystem(lmp)

    # get system & atomic properties
    current_system.get_system_properties()
    current_system.get_atomic_properties()
    current_system.compute_mass_per_type()

    if args.forcefield is not None:

        try:
            ff_atom_type, ff_mass_fraction, ff_partial_charge = np.loadtxt(args.forcefield, \
                    unpack=True, skiprows=1, dtype=[('f0', int), ('f1', float), ('f2', float)])
        except Exception:
            print(f"There was an error while parsing the force field file. \
                    The exception raised is {Exception}. The script will exit!")
            sys.exit()
        if len(ff_atom_type) != n_new_types:
            print("The number of atom types in the force field file differs \
                    from the number of shell types")
            sys.exit()

        polarize(current_system, args, ff_atom_type, ff_mass_fraction, ff_partial_charge)

    else:

        # check that the system has bonds
        if current_system.nbonds == 0:
            print("There are no bonds in the system. The script will exit")
            sys.exit()

        current_system.depolarize(args.shell_atoms)

    # write the output file
    lmp.command(f"write_data {args.output_file} nocoeff ")

    # close the lammps instance
    lmp.close()

if __name__ == '__main__':
    core_shell()
