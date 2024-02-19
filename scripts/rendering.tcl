
# Set the reference structure
set reference_structure "/Users/giuseppegambini/Desktop/Multiscale_Recoverin/Martini_elnet/dynamic/AA_sim/ref_nowat.pdb"

# Load the reference structure
mol new $reference_structure

# Represent the protein
mol representation NewCartoon
mol color ColorID 0
mol selection "protein"
mol material Opaque
mol addrep top

# Calcium ion representation
mol representation VDW
mol color ColorID 4
mol selection "resid 301 or resid 302 and name CAL"
mol material Opaque
mol addrep top

# Set up visualization settings
display depthcue off
display rendermode GLSL
display ambientocclusion on
display antialias on
display projection Orthographic
display nearclip set 0.1
display update

# Set the background color (optional)
color Display Background white

# Update the display to apply the changes
display update

# Plot box
# pbc box

# Display the visualization
display resetview



    