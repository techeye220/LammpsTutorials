# Simulation of a DOPC membrane with the ELBA potential

This tutorial describes how to simulate and analyse a simulation of a DOPC membrane described with the ELBA potential using the [LAMMPS software](http://lammps.sandia.gov/). 

It is intended for new users and such I will go into depth on how LAMMPS works. To this end, I will make extensive references to the LAMMPS [Manual](http://lammps.sandia.gov/doc/Manual.html).

It is assumed that you have downloaded the `Scripts` repository from my Github account to a location referred to as `$SCRIPTS`. It is also assumed that you have installed my copy of LAMMPS to a location referred to as `$LMPDIR`. Both of these things are described in the parent README file.

### LAMMPS input

As we now have a version LAMMPS it is time to setup a simulation of DOPC. We will start from an atomistic representation of the membrane that we will convert to the coarse-grained ELBA model. The file `step5_assembly.pdb` contains the atomistic membrane as assembled by the [CHARMM Membrane builder](http://www.charmm-gui.org/?doc=input/membrane) and the file `step5_assembly.str` contains box information (among other things).

The PDB file can be visualized with for instance VMD

    vmd step5_assembly.pdb

it contains 128 DOPC lipids and 40 water molecules per lipid. 

First, we will convert this to a CG representation with the `aa2cg.py` script. This script needs information about the size of the simulation box that we can find in the `step5_assembly.str`. The lines

     SET A        = 66.7892
     SET B        = 66.7892
     SET C        = 80.96

holds the box size in the x, y and z dimensions. Now type

    python2.7 $SCRIPTS/Lammps/aa2cg.py step5_assembly.pdb -b 66.7892 66.7892 80.96 -i forcefield.elba -o dopc_elba

The `-b` flag sets the box information, the `-i` flag uses a file that we will discuss more below and the `-o` flag sets a label for the output files. The script will create a PDB file with the CG system, `dopc_elba.pdb` that can be visualized with VMD

    vmd dopc_elba.pdb

#### LAMMPS data file

The script will also create a LAMMPS *datafile*, `data.dopc_elba`, that you can read more about [here](http://lammps.sandia.gov/doc/read_data.html). It contains information that LAMMPS needs to run a simulation and can be very complicated. In this tutorial it contains box information, atom definitions and topology information (bond and angles). 

The head of the file should looke something like this

    7040 atoms
    1792 bonds
    2048 angles
    0 dihedrals
    0 impropers

    6 atom types
    5 bond types
    5 angle types

and tell LAMMPS how many atoms, bond and dihedrals there are in the system as well as the number of atom, bond and angle types. The next lines

    -33.39460        33.39460 xlo xhi
    -33.39460        33.39460 ylo yhi
    -40.48000        40.48000 zlo zhi

contains the box information as we specified on the command line above.

Then follows the definition of the atoms in the system.

    Atoms
    
         1  2      6.58200     -0.39600     19.68400    2  0.70000      0.00000      0.00000      0.00000      5.40000      1.10000
         2  3      2.85200      1.64900     18.23000    2 -0.70000      0.00000      0.00000      0.00000      5.20000      1.20000 
    ...

these lines can be very complicated and depends on what type of "atom" you have in your system. You can read more about it [here](http://lammps.sandia.gov/doc/atom_style.html) and [here](http://lammps.sandia.gov/doc/read_data.html). For this simulation the columns have the following meaning

Column(s) | Data
----------| ----
1 | Serial number
2 | Atom type
3-5 | Cartesian (x,y,z) coordinates
6 | Molecule ID number
7 | Charge
8-10 | Dipole moment
11 | Diameter
12 | Density

After the atom definition, you will find the bond definitions

    Bonds
    
         1  1      1      2 
         2  2      2      3 
    ...

the first column is just a serial number, the second column is the bond type and finally, the third and fourth columns are the serial number of the atoms that are bonded together.

Similarly, the datafile also contain angle definitions

    Angles
    
         1  1      1      2      3 
         2  2      2      3      4 
    ...

the format follows the one for the bond definitions and the third to fith columns are the serial number of the atoms forming the angle.

#### LAMMPS include file

You saw in the beginning of the datafile that we defined 6 atom types, 5 bond types and 6 angles types and when we defined atoms, bonds and angles we referred to these types. The force field values of these parameters could be defined in the datafile, but because they are a part of the force field definition and thus will not change from simulation to simulation, we fill put them in an *include* file. This file, `forcefield.elba` will then be "included" or read by LAMMPS to define the necessary force field parameters for the simulation.

The first lines of this file

    pair_style lj/sf/dipole/sf 12.0
    special_bonds lj/coul 0.0 1.0 1.0
    bond_style harmonic
    angle_style hybrid cosine/squared dipole

set the "style" of the pair, bond and angle potential, i.e. what mathematical form they will take. The pair style will be a shifted-force Lennard-Jones for the van der Waals interaction and a shifted-force dipole/charge potential for the electrostatics. The final `12.0` sets the non-bonded cut-off to 12 A. Furthermore, the `0.0 1.0 1.0` numbers in the `special_bonds` command indicate scaling of 1-2, 1-3 and 1-4 interactions. You can read more about pair styles [here](http://lammps.sandia.gov/doc/pair_style.html)

The bond style will be a simple harmonic function. 

Finally, the angle style will be a hybrid, i.e. one or more possible functional forms. This is an example of the very flexible nature of LAMMPS were an angle potential can be of different styles. The `cosine/squared` will be used for the standard valance angles whereas the `dipole` style will be used to restrain dipoles in the lipids.

The rest of the `forcefield.elba` file sets 1) the mass, diameter and Lennard-Jones parameters for the different atom types, 2) set the coefficients for the bond types and 3) the coefficients for the angles types.


#### LAMMPS input file

The final piece of input is the actual input file that LAMMPS reads to perform a simulation (or a set of simulations). You can read more about the input file [here](http://lammps.sandia.gov/doc/Section_commands.html). The important thing to notice here is that LAMMPS read this file one line at a time and execute this line as a command. Therefore, the order of commands are very important.

Let's go through the input file in detail

The first section

    units	real
    atom_style	hybrid angle dipole sphere 
    read_data 	data.dopc_elba
    include 	forcefield.elba
    velocity	all create 0.0 87287 

setup up the simulation. It sets the units and the style of the atoms before it reads our two other pieces of input, the datafile and the include file. Then it will initialise the velocities for all particles.

The second section defines a few variables

    variable	nLips equal 128 # total number of lipids
    variable	nWats equal 5120 # total number of waters 
    variable	watVol equal 30.0 # water molecular volume (~30 A^3)

LAMMPS has the capabilities to create a wide range of variable that can be quite complicated. You can read more about it [here](http://lammps.sandia.gov/doc/variable.html). Here we define some constants so that we later can calculate the area and volume per lipid. 

The third section defines some groups

    group		lip type 2 3 4 5 6 
    group		head type 2 3 
    group		wat type 1
    group		chol type 2
    ...

groups are useful to define so that we can perform various operations on them. This will be clearler below. You can read more about groups [here](http://lammps.sandia.gov/doc/group.html).

Furthermore, we will set the time step to 10 fs with

    timestep	10

The first piece of action takes place with the command

    minimize        1.0E-4 1.0E-6 100 1000

that tell LAMMPS to minimise the system for 100 steps. You can read more about the command [here](http://lammps.sandia.gov/doc/minimize.html)

After that we will run a short simulation in the NVT ensemble. Therefore, we need to set and integrator and a thermostat. This is accomplished in LAMMPS through so-called [*fixes*](http://lammps.sandia.gov/doc/fix.html) that are actions performed on a selection of atoms every timestep. The general formula for a fix command is

    fix *name* *atom_selection* *name_of_fix* *settings*

Therefore, the lines

    fix   integrate all nve/sphere update dipole
    fix   thermo all langevin 303 303 1000 9 omega yes zero yes

first tell LAMMPS to use the fix `nve/sphere` on the selection `all` that has the effect to propagate the motion of all particles according to velocity Verlet integrator. Second, it tells LAMMPS to apply the fix `langevin` to all atoms that has the effect that the system will be simulated at 303 K using a Langevin thermostat.

To run the simulation we simply use the command

    run		5000


After the NVT ensemble, we want to run a short equilibration in the NPT ensemble. Therefore, we will introduce a fix that controls the pressure in the simulation

    fix		baro all press/berendsen aniso 1 1 1000 couple xy modulus 21740
    run		50000

As you see it is very easy to run several simulations after each other in LAMMPS without the need to create several input files.

Finally, we will run a 10 ns simulation where we will output some interesting data for analysis. 

Some output will be realised through a combination of variables and fixes and a trajectory will be written out with the *dump* command. These commands are rather complicated and can be studied in detail in the LAMMPS manual. In short, we will create a range of files with number density for the different beads and we will output the area and volume per lipid.


### Running and analysing the simulation

To run LAMMPS using 8 processors type

    mpirun -np 8 ${LMPDIR}/src/lmp_soton < in.sim

where `$LMPDIR` is the installation directory of LAMMPS. It will take a couple of hours to complete.

It will produce a number of output files. `apl.dat` and `vpl.dat` contains the trajectory of the area and volume per lipid. You can average them with a simple script that is provided

    python av.py < apl.dat
    python av.py < vpl.dat

(this script assumes you have Numpy installed).

Furthermore, the simulation will produce a number of files starting with `numDens`. These are number densities for the different bead types along the membrane normal. You can copy them into an Excel sheet an plot them to see the density profile of the membrane.

Finally, you can visualise the entire trajectory with VMD using

    vmd dopc_elba.pdb sim.dcd

as you will notice the lipids are broken over the central simulation box and thus you will see very long bonds drawn in VMD. It is therefore best to use a VDW representation of the beads.

This tutorial has now showed you how to setup a simple membrane simulation LAMMPS. It has showed you what the different pieces of input are and how to prepare and understand them. Finally, it has showed you some simple and straightforward analysis.


