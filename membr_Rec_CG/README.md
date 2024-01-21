
# First attempt

The objective is to generate the CG protein with `martinize2` and then use `insane` to insert the myristoyl group in the membrane. 

I can use the gromacs option `-membed` in `mdrun` to insert the protein inside the membrane. The method works by first artificially shrinking the protein in the xy-plane, then it removes lipids that overlap with that much smaller core. 

The outline of the protocol is the following:
1. Martinize2 to coarse-grain the protein structure (specify the parameters for the myristoyl group)
2. Generate a membrane system using `insane`
3. Use again `insane` to place the myristoyl inside the membrane
4. Add water molecules to solvate the system and add ions to neutralize
5. Perform an energy minimization to relax the system and remove any steric clashes.



# Point 1: martinize the protein
First of all, we isolate the protein in `Rec.pdb`. I can do this simply with VMD (or also trjconv). 

Then we have to extract the secondary structure, which has to be fed to Martinize2. To do this, we use the online tool http://bioinformatica.isa.cnr.it/SUSAN/DSSP-web/, copy the `.pdb` file and then wait for the output sequence. We proceed with We proceed with 
```
martinize2 -f Rec.pdb -o Rec.top -x Rec_cg.pdb -ss ccccccSccSHHHHcTTccccSSHHHHHHHHHHHHTTSGGGcEEHHHHHHHHHHHScSSccHHHHHHHHHHHcSScSSEEcHHHHHHHHHHHSccSSSccHHHHHHHHcSScSSEEcHHHHHHHHHHHHHHScHHHHHHScSTTTSHHHHHHHHHHHTTccTTccEEHHHHHHHHHHcHHHHHHTcccHHHHHHHHHHccc -p backbone -ff-dir charmm36-feb2021_vm.ff -ff martini3001 -nt -map-dir ../parameters/parameters/myristoyl
```
The specification `-nt` is important because it makes the termini neutral. 

# Point 2: generate the membrane 
Insane is a CG building tool that generates membranes by distributing lipids over a grid

# Point 3: place the protein inside the membrane

# Point 4: solvation

# Point 5: energy minimization

