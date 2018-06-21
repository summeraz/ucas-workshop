import MDAnalysis as mda

output = mda.Universe('output.xyz', format='XYZ')
output.atoms.write('output.pdb')
