import pdb_objects as pdbo
parser = pdbo.PDBParser()
parser.parse('ADHE_____.pdb')
pdb = parser.parse('ADHE_____.pdb')
model = pdb.model0
dihs = model.get_dihedrals()

for chain, ds in dihs.items():
	f = open('dihedrals_chain_{}.csv'.format(chain),'w')
	f.write('{}\t{}\t{}\t{}\n'.format('RESNUM', 'RESNAME', 'PHI', 'PSI'))
	for (num, name), (phi, psi) in ds:
		f.write('{}\t{}\t{}\t{}\n'.format(num, name, phi, psi))
	print 'Dihedrals of chain "{}"" written to "dihedrals_chain_{}.csv"'.format(chain,chain)
	f.close()