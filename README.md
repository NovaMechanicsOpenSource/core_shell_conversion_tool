
**L-CS-RI: A LAMMPS core-shell to/from rigid-ion configuration converter**

**Introduction**

The _core-shell model_ is an easy way of introducing polarizability into classical models, with applications primarily in crystalline and amorphous inorganic solids. In this model, an atom is split to two particles, one forming the core and a second forming the shell which corresponds to the relative motion of electrons to atom's nucleus. Core-shell models do not require that all particles in a system are described by a core-shell geometry.

LAMMPS can carry out simulations based on the core-shell model using a variety of interatomic potentials. These potentials are included in the CORESHELL package; two examples are the suitably modified Buckingham and 12-6 Lennard-Jones potentials. Nevertheless, LAMMPS cannot convert data- or restart files of rigid-ion configurations to corresponding core-shell ones and vice versa. Thus, users need to rely on external tools or develop their own custom codes. The **L-CS-RI** tool offers such a converter based on the LAMMPS library interface. It comes with a command-line interface (CLI), making it suitable for use in high throughput simulation workflows.

The steps for converting a rigid-ion to a core-shell configuration are: (i) splitting the mass of each atom type based on the force field in the mass of the core and the shell types, (ii) splitting the partial charge of each atom by adding the charge indicated by the force field, (iii) creating bonds between the core and the shell particles, and (iv) creating the coordinates of the shell atoms by randomly displacing the coordinates of the rigid ion by a small amount.

The steps for converting a core-shell to a rigid-ion configuration are: (i) removing the shell particles, (ii) deleting all bonds that the shell particles are involved, (iii) adding the mass of the shell particles to the core particles where they are bonded, and (iv) adding the partial charge of the shell particles to their core particles.

**Description of the code**

The **L-CS-RI** tool consists of two files: the particle\_system.py, which contains the ParticleSystem class, and the core\_shell.py. A description of the class and its methods is shown in Tables 1 and 2, respectively. **L-CS-RI** requires NumPy 1.13.3 [1], ctypes 1.1.0, which is part of the Python standard library, and LAMMPS shared library version 22 Dec 2022 or later [2]. **L-CS-RI** has been tested with python 3.6.9.

A rigid-ion configuration can be converted to a core-shell one by issuing the command:

_python core\_shell.py__input\_file output\_file list\_of\_particle\_types --restart\_file --atom\_style {ATOM\_STYLE} --extra\_keywords {EXTRA\_KEYWORDS} --forcefield {FORCEFIELD}_

The inverse procedure of merging core-shell particles to rigid-ions can be accomplished by:

_python core\_shell.py__input\_file output\_file list\_of\_particle\_types --restart\_file --atom\_style {ATOM\_STYLE} --extra\_keywords {EXTRA\_KEYWORDS}_

The command line interface of the **L-CS-RI** tool has three positional and four optional arguments:

1. _input\_file_: A string with the absolute path and filename of the input file that will be processed. The file should be a LAMMPS data file containing the initial configuration of the system of interest. The code ignores any information regarding the functional form and the parameters of the force field that might be present in the file.
  1. If the input file is a restart file, i.e. binary, then the --_restart\_file_ flag should be used.
  2. If specific keywords are required when reading the file, they can be provided as string using the _--extra\_keywords_ flag.
2. _output\_file_: A string the absolute path and filename of the output file with the generated configuration. The output file is always a text file complying with the format of LAMMPS data files.
3. _list\_of\_particle\_types_: a list of integers corresponding to the particle types that are either (i) rigid ions and will be split into core and shell particles or (ii) shell particles that will be merged with their core particles to form the rigid ions.

The LAMMPS _atom\_style_ attribute is by default "full". If the user wishes to specify a different style, which supports partial charges and bonds per particle type, then the _--atom\_style_ flag followed by the desired _style_ should be specified.

If the configuration will be polarized, i.e. shell particles will be added to the system, then the _--forcefield_ flag followed by a string indicating the full name and path to the force field file must be provided. The format of the force field file is:

- The first line is a considered to be a comment line and it is ignored.
- The rest of the lines should contain three columns. The first column is the atom type (an integer number) corresponding to the rigid-ion that will be split and become a core particle. The second is the mass fraction that will be retained in the core particle (real number between 0.0 and 1.0). The third is the partial charge (real number) that will be added to the core particle and will be subtracted from the shell particle.

The shorthand notation of the optional arguments is (i) _-a_ for _--atom\_style_, (ii) _-e_ for _--extra\_keywords_, (iii) _-f_ for _–forcefield_ and _-r_ for _--restart\_file_.

A short description of the tool and its usage can be obtained by issuing the command:

_python core\_shell.py --help_

**Examples**

The first example is the conversion of a NaCl core-shell configuration to a rigid-ion one. The force field has been derived by Mitchell and Fincham [3]. The input file is already in the LAMMPS distribution under the examples/coreshell directory. The command for the conversion is

_python core\_shell.py data.coreshell data.rigid\_ion 3 4 -e "fix csinfo NULL CS-Info"_

The second example is the conversion of the previously obtained core-shell configuration to the initial rigid-ion one:

_python core\_shell.py data.rigid\_ion data.new\_coreshell 1 2 -f core\_shell.ff_

The content of the _core\_shell.ff_ force field file is as follows:

_# force field file
 1 0.9 -0.5005
 2 0.9 -2.5056_

**References**

1. Harris, C.R., Millman, K.J., van der Walt, S.J. et al. "Array programming with NumPy". Nature 585 (2020) 357.
2. A. P. Thompson, H. M. Aktulga, R. Berger, D. S. Bolintineanu, W. M. Brown, P. S. Crozier, P. J. in 't Veld, A. Kohlmeyer, S. G. Moore, T. D. Nguyen, R. Shan, M. J. Stevens, J. Tranchida, C. Trott, S. J. Plimpton "LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales" Comp Phys Comm 271 (2022) 10817
3. P J Mitchell and D Fincham "Shell model simulations by adiabatic dynamics" J. Phys.: Condens. Matter 5 (1993) 1031

Table 1: Description of the ParticleSystem class in the particle\_system.py file

| Attributes of the ParticleSystem class | Short description |
| --- | --- |
| lmp\_pointer | Instance of the LAMMPS Python class |
| Ntypes | Number of atom types |
| Nbondtypes | Number of bond types |
| natoms | Number of atoms |
| nbonds | Number of bonds |
| mass\_per\_type | Vector of ntypes length containing the mass of each atom type |
| atom\_id | Vector of natoms length containing the id of each particle |
| atom\_type | Vector of natoms length containing the type of each particle |
| atom\_mass | Vector of natoms length containing the mass of each particle |
| atom\_partial\_charge | Vector of natoms length containing the partial charge of each particle |
| atom\_coordinates | Vector of 3\*natoms length containing the coordinates of each particle |
| bonds | Vector of 3\*nbonds length containing the bond type and particles of each bond |
| Methods of the ParticleSystem class |
| \_\_init\_\_ | initialize the object variables |
| get\_system\_properties | get system information from lammps object |
| get\_atomic\_properties | get atomic information from lammps object |
| compute\_mass\_per\_type | compute the mass of each atom type from lammps object |
| set\_atomic\_properties | send the information to lammps for several properies |
| depolarize | remove the shell particles |

Table 2: Functions in the core\_shell.py file

| core\_shell | Main function for adding/removing shell particles from a lammps data/restart file |
| --- | --- |
| create\_cli | Create an elementary CLI based on the argparse module & perform initial consistency checks |
| polarize | Create the shell particles |
