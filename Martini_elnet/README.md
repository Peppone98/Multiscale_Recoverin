
# Martini + elastic network 

We add a supplementary level of constraints on the protein, namely, harmonic springs between the beads. More precisely, we add to the standard Martini topology extra harmonic bonds between non-bonded beads based on a distance cut-off. 

Martinize2 does not generate elastic bonds between $i\rightarrow i+1$ and $i\rightarrow i+2$ backbone beads, as those are already connected by bonds and angles.

Martinize2 will generate harmonic bonds between backbone beads if the option `-elastic` is set. The bond strengths independent of the bond length
```
martinize2 -f nmRec.pdb -o nmRec_el.top -x nmRec_cg_el.pdb -ss ccTTcScccSSTTTSGGGcScScHHHHHHHHHHHHHHSGGGcBcHHHHHHHHHHHSTTcccHHHHHHHHHHHcSSSSSSBcHHHHHHHHHHHccSSSSccHHHHHHHHcSScSScBcHHHHHHHHHHHHHHScHHHHHTSSTTccSHHHHHHHHHHHTTccTTccBcHHHHHHHHHHcHHHHHHHcccHHHHHHHHHcccc -p backbone -ff martini3001 -nt -elastic -ef 500.0 -el 0.5 -eu 0.9
```
Here are some reasons why the cutoff distance in an elastic network is often set to a value around 0.9 nm or smaller: 

Now I manually edit the `.top` file in order to include the two structural Calcium ions. From now on, look at the `system_el.top` file. To avoid unwanted errors of numeration, I labelled the two `CA` residues as residues 202 and 203. The CG coordinates are saved in `system_cg_el.pdb`.

I create the file `system_cg_el.pdb` as a copy of `nmRec_cg_el.pdb` and I add the two ion coordinates to this.

Importantly, if you want to rename the ugly `molecule_0.itp`, then you must also modify the file itself after `[ moleculetype ]`. Then of course you have to modify the topology file accordingly. 

# Short minimization in vacuum 
First, we define the simulation box. As suggested by Alessio, I enlarge a little bit the box:
```
gmx editconf -f system_cg_el.pdb -box 12 12 12 -o system_cg_el.gro
```
In the `minimization.mdp` file we leave 10 kJ/mol (the default value for `emtol` in gromacs). Remember to correct the `system.top` with the name of the `martini.itp` you are using. 
```
gmx grompp -p system_el.top -f mdp/minimization.mdp -c system_cg_el.gro -o tpr/min_vacuum.tpr -r system_cg_el.gro
```
```
cd min_vacuum
gmx mdrun -s ../tpr/min_vacuum.tpr -v 
```
# Solvate the system
Gromacs needs to know the van der Waals radii of each atom, namely, the radius of an imaginary hard sphere representing the distance of closest approach for another atom. Think at the van der Waals equation of state: there you have terms accounting for the efective volume occupied by the atoms. For instance, the $r_{\text{vdw}}$ of Carbon is 0.17 nm. Since we are introducing new atoms with respect to the ones listed in the default file `vdwradii.dat` of gromacs, we have to specify a global vdw radius. 
```
gmx solvate -cp min_vacuum/confout.gro -cs water.gro -radius 0.21 -o system_solv.gro
```
The value 0.21 typically avoids clashes because gromacs will then fill the space available to water by trying to be gentle with the distances between beads (of both water and protein). In fact, if two beads are generated too close, this situation may correspond to an overlap in the AA system, which is clearly unwanted.

Gromacs will return you a warning because you are not passing a `.tpr` file (so he doesn't know how to assign the vdw radius).

# Adding the ions
You first need to update the topology to reflect the added water (`W        13336`). And then you must also add the topology file of the solvent (all these operations aren't done automatically as in groamcs).
```
gmx grompp -f mdp/minimization.mdp -c system_solv.gro -p system_el.top -o tpr/ions.tpr
```
```
gmx genion -s tpr/ions.tpr -o solv_ions.gro -p system_el.top -pname NA -nname CL -conc 0.15 -neutral
```
The potassium ion is not present in the Martini force field for ions, so I use sodium instead. Long range electrostatic interactions are absent, and for small ions the first hydration shell is considered an implicit part of the CG ion.

## Minimization in solution
The only difference is that here we do not have the protein in vacuum anymore, but the whole system. 
```
gmx grompp -p system_el.top -c solv_ions.gro -f mdp/minimization.mdp -o tpr/min_solution.tpr -r solv_ions.gro
```
```
cd min_solution
gmx mdrun -s ../tpr/min_solution.tpr -v
```

# Equilibration run: NPT
Most simulations are numerically stable with dt=40 fs. Remember to include the line `define = -DPOSRES` in the `.mdp` file. The Martinize2 script already inserted the restraints on the backbone atoms in `nmRec.itp`. You may also need to use `refcoord_scaling = com`. 
```
gmx grompp -p system_el.top -c min_solution/confout.gro -f mdp/equilibration.mdp -o tpr/equilibration.tpr -r min_solution/confout.gro
```
```
cd equilibration
gmx mdrun -s ../tpr/equilibration.tpr -v
```

# Production run
```
gmx grompp -p system_el.top -c equilibration/confout.gro -f mdp/dynamic.mdp -o tpr/dynamic.tpr 
```


# NMR data 
Using NMR data to determine the spring constant for an elastic network in a coarse-grained simulation can be a valuable approach to better represent the dynamics of a protein. 

Perform short simulations (just a few nanoseconds) with different spring constants and compare the simulated dynamics with the NMR-derived dynamics. 


