# LAMMPS Tutorials

This repositoty contains LAMMPS tutorials for simulations with the ELBA model.

For now it contains the following tutorials

1.) Membrane - simulation of a DOPC membrane with the ELBA potential
2.) ElbaSolv - solvation free energy calculation of an ELBA water in a box of ELBA water molecules



It is assumed that you have downloaded the `Scripts` repository from my Github account. This can be accomplished with

    git clone https://github.com/SGenheden/Scripts

In the tutorials, I will refer to the location of the `Scripts` folder with `$SCRIPTS`.

### Installing LAMMPS

We will start with installing LAMMPS. We will use my own version of LAMMPS hosted on Github. It contains some of my modifications to the code and it will be required by some ELBA tutorials. Some of the modifications will eventually be released in the official version of Lammps.

Clone my lammps repository from Git hub

    git clone https://github.com/SGenheden/lammps

Next, we will install the colvars library that is needed for some simulations

    cd lammps/lib/colvars
    make -f Makefile.g++

Now we are ready to install LAMMPS. By default LAMMPS will only contain a limited set of functional. The rest of the functionality can be accessed by installing a set of packages. You can read more about packages [here](http://lammps.sandia.gov/doc/Section_start.html#start_3). Therefore, we will start to add the necessary packages for ELBA simulations

    cd ../../src
    make yes-dipole yes-kspace yes-manybody yes-misc yes-molecule yes-rigid yes-user-colvars
    cp USER-MISC/angle_dipole.* USER-MISC/pair_lj_charmm_coul_long_14.* USER-MISC/pair_lj_sf_dipole_sf.* .

The last command add a few bits from the USER-MISC package. This is a large package with user-contributed code and it is unnecessary to install all of it.

To compile LAMMPS we need an appropriate make file. One that work with the OpenMPI and Intell compilers on Southampton machines/clusters can be found in this repository. Copy the file `Makefile.soton` to the `lammps/src/MAKE` folder and then execute

    make -j 8 soton

The installation should end with something like this

    size ../lmp_soton
       text        data     bss     dec     hex filename
    8478121      176400   18176 8672697  8455b9 ../lmp_soton
