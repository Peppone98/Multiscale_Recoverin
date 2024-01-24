
# Short comparison with the all atom (AA) trajectory

1. Distance between the six beads of the coordinant oxygens and the calcium ion
2. $\Delta$RMSD between AA and CG 
3. $\Delta$RMSF between AA and CG 


The first thing I have to do is to provide the correct alignment of the protein. I want to center the system on the protein and the calcium ion. First I have to define an index 
```
gmx make_ndx -f ../dynamic/equilibration.pdb
```
Then I center the trajectory 
```
gmx trjconv -s ../tpr/dynamic.tpr -f ../dynamic/traj_comp.xtc -pbc nojump -center -o traj_centered.xtc -n index.ndx
```
I did the fit with respect to the protein CG sites. 

```
gmx rmsdist -f ../dynamic/traj_comp.xtc -s ../tpr/dynamic.tpr -o RMSD.xvg
```

```
gmx rmsf -f ../dynamic/traj_comp.xtc -s ../tpr/dynamic.tpr -o RMSF.xvg
```



