'''This modules defines the class for particle systems '''

import numpy as np

class ParticleSystem:
    '''Class to store data of a lammps particle system'''

    def __init__(self, lmp_object):
        '''initialize the object variables'''
        self.lmp_pointer = lmp_object
        self.ntypes = self.lmp_pointer.extract_setting("ntypes")
        self.nbondtypes = None
        self.natoms = None
        self.nbonds = None
        self.mass_per_type = [0] * self.ntypes
        self.atom_id = None
        self.atom_type = None
        self.atom_mass = None
        self.atom_partial_charge = None
        self.atom_mol = None
        self.atom_x_coordinates = None
        self.atom_y_coordinates = None
        self.atom_z_coordinates = None
        self.bonds = None

    def get_system_properties(self):
        '''get system information from lammps object'''
        self.nbondtypes = self.lmp_pointer.extract_setting("nbondtypes")
        self.natoms = self.lmp_pointer.get_natoms()
        self.nbonds = self.lmp_pointer.extract_global("nbonds")

    def get_atomic_properties(self):
        '''get atomic information from lammps object'''
        self.bonds = list(self.lmp_pointer.gather_bonds()[1])

        self.lmp_pointer.command("variable atom_id atom id")
        self.atom_id = np.array(list(self.lmp_pointer.extract_variable("atom_id")), dtype=int)

        self.lmp_pointer.command("variable atom_type atom type")
        self.atom_type = np.array(list(self.lmp_pointer.extract_variable("atom_type")), dtype=int)

        self.lmp_pointer.command("variable charge atom q")
        self.atom_partial_charge = np.array(list(self.lmp_pointer.extract_variable("charge")))

        try:
            self.lmp_pointer.command("variable molecule atom mol")
            self.atom_mol = np.array(list(self.lmp_pointer.extract_variable("molecule")), dtype=int)
        except Exception:
            print("It was not possible to retrieve the molecule id from LAMMPS")
            self.atom_mol = None

        self.lmp_pointer.command("variable atom_mass atom mass")
        self.atom_mass = np.array(list(self.lmp_pointer.extract_variable("atom_mass")))

        self.lmp_pointer.command("variable atom_x atom x")
        self.atom_x_coordinates = np.array(list(self.lmp_pointer.extract_variable("atom_x")))

        self.lmp_pointer.command("variable atom_y atom y")
        self.atom_y_coordinates = np.array(list(self.lmp_pointer.extract_variable("atom_y")))

        self.lmp_pointer.command("variable atom_z atom z")
        self.atom_z_coordinates = np.array(list(self.lmp_pointer.extract_variable("atom_z")))

    def compute_mass_per_type(self):
        '''compute the mass of each atom type from lammps object'''
        for itype in range(0, self.ntypes):
            if np.count_nonzero(self.atom_type == (itype+1)) > 0:
                self.mass_per_type[itype] = np.average(self.atom_mass[self.atom_type == (itype+1)])
            else:
                self.mass_per_type[itype] = np.finfo(float).eps # a small number

    def set_mass_per_type(self):
        '''send the information to lammps for mass'''
        for itype, imass in enumerate(self.mass_per_type):
            self.lmp_pointer.command(f"mass {itype+1} {imass}")

    def depolarize(self, shell_types):
        '''remove the shell particles'''

        # fix the partial charges & the masses of the remaining particles
        for iat, jat in zip(self.bonds[1::3], self.bonds[2::3]):
            idx = np.where(self.atom_id==iat)[0][0]
            idy = np.where(self.atom_id==jat)[0][0]
            if self.atom_type[idx] in shell_types:
                self.atom_partial_charge[idy] += self.atom_partial_charge[idx]
                self.atom_partial_charge[idx] = 0
                self.atom_mass[idy] += self.atom_mass[idx]
                self.atom_mass[idx] = np.finfo(float).eps # a small number
            elif self.atom_type[idy] in shell_types:
                self.atom_partial_charge[idx] += self.atom_partial_charge[idy]
                self.atom_partial_charge[idy] = 0
                self.atom_mass[idx] += self.atom_mass[idy]
                self.atom_mass[idy] = np.finfo(float).eps  # a small number

        for iat, new_q in zip(self.atom_id, self.atom_partial_charge):
            self.lmp_pointer.command(f"set atom {iat} charge {new_q}")

        # set the proper values for the mass
        self.compute_mass_per_type()
        self.set_mass_per_type()

        # remove the shell particles  & reset the atom ids
        shell_string_list = [str(x) for x in shell_types]
        self.lmp_pointer.command(f"group shell_particles type {' '.join(shell_string_list)}")
        self.lmp_pointer.command("delete_atoms group shell_particles bond yes")
        self.lmp_pointer.command("reset_atoms id sort yes")
