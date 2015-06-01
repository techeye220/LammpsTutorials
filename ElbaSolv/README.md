# Solvation free energy of water with the ELBA potential

This tutorial describes how to setup, run and analyse the solvation free energy simulation using Lammps and the coarse-grained ELBA potential.

It assumes that you have cloned and installed the Lammps version in my GitHub repository.

### Creating a box of water

As a start, we will create a box of atomistic waters that surrounds a single, central water molecule. This central water molecule will then be the subject of the solvation free energy calculation. 

You can skip this step if you do not have my python scripts, as the CG box is provided with this tutorial.

We will create this box with *tleap* from the AmberTools

    cat << EOF > leapcom
    source leaprc.ff99SB
    x=copy TP3
    solvatebox x TIP3PBOX 10
    savepdb x box.pdb
    quit
    EOF

    tleap -f leapcom
    
and then we can create it to a CG representation with the *aa2cg.py* script. Write down the box dimensions that *tleap* writes out.

    sed -i "s/O  /OH2/" box.pdb
    python2.7 $SCRIPTS/Lammps/aa2cg.py -f box.pdb -o box_elba -i forcefield.elba -b 26.782 27.009 26.192
    
this will create *box_elba.pdb* that you can use to visualize the system and *data.box_elba* that is the configuration+topology input to lammps. 

The force field terms are provided with this tutorial in *forcefield.elba*

### Exploring and editing the datafile and forcefield file

**The datafile**

The datafile contains configuration and topology (bonds, angles etc). As we have a CG system of only water beads, there is no topology. Each atom line in the datafile describes one atom and its property. 

The information in each column depends on the *atom_style*. We will use a *hybrid molecular dipole sphere* style, which means that the column contains the following information: serial number, atom type, Cartesian coordinates, molecule id, charge, dipole vector, diameter and density. 

You will see that each atom has an atom type 1. However, we need to be able to distinguish between the central atom, which is our solute and the other water molecules. Therefore, change the atom type of the first atom to 2. The line should read something like this

    1  2     -0.18400     -0.03000     -0.12200    1 ...

Also change the line 

    1 atom types
    
to
  
    2 atom types   

so that Lammps understand to look for two atom types.     

**The force field file**

The force field file is an inclusion file, a snippet of a Lammps input file that is put in a separate file for convenience. At the moment it only contains a *mass* and a *pair_coeff* command

We need to add mass and pair coefficients two our second, artificial atom type by adding the following lines

    mass    2   40.000  # wat
    
and 

    pair_coeff    1    2      0.550  3.050 # wat-wat
    pair_coeff    2    2      0.550  3.050 # wat-wat
    
so that Lammps know the mass of atom type 2 and how atom type 2 and atom type 1 interacts with each other.

### Equilibration

It is good to start with a short equilibration. Here will run a 50 ps simulation in the NPT ensemble, using a 10 fs timestep. 

The input script *in.equil* is commented fairly well and should be easy to understand. We will keep the central atom (our solute) fixed at the centre of the box.

To run the simulation type something like this
    
    mpirun -np 8 $LMPPATH/src/lmp_soton < in.equil > out.equil
    
where `$LMPPATH` is the installation directory of Lammps. The simulation should only take a few minutes. You can visualize the simulation with for instance VMD.

    vmd box_elba.pdb equil.dcd

The equilibration will produce new datafile that will be the starting point for the solvation free energy calculation.
    
### Solvation free energy

Now we have an equilibrated box of ELBA waters and we are ready to run the solvation free energy calculation. 

The procedure very much follows the one in Orsi, M., W. Ding, M. Palaiokostas. *J. Chem Theory Comput*, **2014**, *10*, 4684.

We will run one single, long simulation and occasionally change *lambda* our scaling parameter. For each value of *lambda*, we will run 25 ps equilibration and 350 ps sampling, and we will sample *lambda* from 0 to 0.96. And scale our system with the fourth power of *1 - lambda*. This is necessary as we cannot use soft-core potentials. We will then make a linear interpolation to *lambda* equal 1.

The input *in.sim* is comment fairly well, but warrants some extra comments. 

The actual scaling takes place with the *fix adapt* command with the lines

    fix adaptSS all adapt 1 pair lj/sf/dipole/sf epsilon 1 2 v_fL &
	                        pair lj/sf/dipole/sf scale   1 2 v_fL scale yes
	                        
which tell Lammps that booth the van der Waals (epsilon) and electrostatics (scale) of the pair style *lj/sf/dipole/sf* should be scaled with the variable *fL*. This variable is initialised with the following three lines of input

    variable stepNow equal step
    variable L       equal (floor((v_stepNow-1)/v_Nwintot)*v_deltaL)
    variable fL      equal (1-v_L)*(1-v_L)*(1-v_L)*(1-v_L)

which defines a variable for the current step, a variable L (*lambda*) which is a function of the current step and the *fL* which is the fourth power of *1 - lambda*. *Nwintot* is the total number of steps for *lambda* value and *deltaL* is the increment in *lambda*, which is set to 0.04. 

The computation of the derivate with respect to *lambda* is carried out with the following lines of input

    compute  PotEngSS  centralwat group/group otherwat pair yes
    variable dPotEngSS equal c_PotEngSS*v_dfL/v_fL
    fix dPotEngSS all ave/time ${Ne} ${Nr} ${Nf} v_L v_dPotEngSS file out.dPotEngSS

the first line computes the interaction between the solute and the other water molecules the second line computes the derivative with respect to *lambda* and the third line prints it to an ouput file *out.dPotEngSS*. The combination of *Ne*, *Nr* and *Nf* ensures that the average is only over the sampling portion and that the output is only once for each *lambda*. The explanation for these variables can be found in the Lammps manual. 

To run the simulation type something like this
    
    mpirun -np 8 $LMPPATH/src/lmp_soton < in.sim > out.sim
    
the simulation will take less than an hour. 

The most important output can be found in *out.dPotEngSS*. This is the derivative of the potential energy with respect to *lambda*. First, we need to find the value of the derivative at *lambda* equals 1 using linear interpolation. Second, we need to integrate this to obtain the free energy of decoupling, i.e. the negative of the free energy of solvation.

A quick-and-dirty python script that uses NumPy is provided

    python integrate.py out.dPotEngSS

should give roughly 6.5 kcal/mol which corresponds well to the experimental hydration free energy of water.
