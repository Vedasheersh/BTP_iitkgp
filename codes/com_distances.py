import pdb_objects as pdb

### This is a very specific code. Better not to use it for general purposes.

def main(ens_pdb_names, res_range_HOX, res_range_flex):
	'''
	Takes care of pdb_parsing, center of mass and distance calculations
	'''

	distances_matrix = []

	## Create parser object
	parser = pdb.PDBParser()

	## Loop over all pdb names present in the ensemble list
	for pdb_name in ens_pdb_names:
		## Create a pdb object (using parser) for the current pdb_name
		pdb_obj = parser.parse(pdb_name)

		## If there are multiple models in the pdb file, there is ambiguity as which model to choose!
		if len(pdb_obj.models) > 1:
			print pdb_obj.models
			model_num = raw_input('PDB file {} has multiple models in it. Which to consider?\
				Enter model number:'.format(pdb_name))
		else:
			## If there is only single model in the pdb, it is numbered 0 by default
			model_num = 0

		## The model object
		_model = pdb_obj.model(model_num)

		## Again, if there are multiple chains in the model choosen, there is ambiguity as which chain! 
		if len(_model.chains) > 1:
			print _model.chains
			chain_name = raw_input('PDB file {} has multiple chains in it. Which to consider?\
				Enter chain name:'.format(pdb_name))
			_chain = _model.chain(chain_name)

		else:
			## If there is only single chain, take it.
			_chain = _model.chains.values()[0]
		
		## Place holder for N, CA, C and O position vectors of all residues of HOX domain
		HOX_heavy_vectors = []
		## Place holder for CA position vectors of all residues of flexible region
		flex_CA_vectors = []

		## Loop through all the residues of chosen chain
		for residue in _chain:
			if res_range_HOX[0]<=residue.serial<=res_range_HOX[1]:
				HOX_heavy_vectors.extend([residue.N.position,residue.CA.position,residue.C.position,residue.O.position])
				# print 'Residue',residue.serial, 'is in HOX.'
			elif res_range_flex[0]<=residue.serial<=res_range_flex[1]:
				flex_CA_vectors.append(residue.CA.position)
				# print 'Residue',residue.serial, 'is in flexible region.'

		## Assuming equal mass for N, C, O ; center of mass position vector is simply the sum of position vectors
		## of position vectors of N,C,O of all the residues.
		ref_vec = HOX_heavy_vectors[0] - HOX_heavy_vectors[0]
		center_of_mass = ref_vec
		for vec in HOX_heavy_vectors[0:]:
			center_of_mass+=vec-ref_vec
		center_of_mass/=len(HOX_heavy_vectors)
		# print center_of_mass, 'is the position vector of center of mass.'

		matrix = [vec.distance(center_of_mass) for vec in flex_CA_vectors]

		distances_matrix.append(matrix)

	return distances_matrix

def plot_results(distances_matrix):
	'''
	Plot the distances from calculated matrix.
	'''
	import matplotlib.pyplot as plt
	matrix_by_res = []
	for res_num in range(len(distances_matrix[0])):
		temp = []	
		for conf_num in range(len(distances_matrix)):
			temp.append(distances_matrix[conf_num][res_num])
		matrix_by_res.append(pdb.np.array(temp))

	avg_matrix = []
	std_matrix = []
	for matrix in matrix_by_res:
		avg_matrix.append(pdb.np.average(matrix))
		std_matrix.append(pdb.np.std(matrix))

	plt.figure()
	plt.errorbar(range(len(avg_matrix)),avg_matrix,yerr=std_matrix,fmt='o')
	plt.plot(avg_matrix,'c',label='average distance')
	plt.plot(std_matrix,'--',label='standard deviation')
	plt.show()

if __name__ == '__main__':

	import argparse
	import json
	parser = argparse.ArgumentParser(description='DISTANCE PLOTS B/W HOX AND EACH RES OF FLEXIBLE REGION')
	parser.add_argument("-hox_range",'--hox_res_range',help="Range of residues in pdb for HOX domain")
	parser.add_argument("-flex_range",'--felx_res_range',help="Range of residues in pdb for flexible region")
	parser.add_argument('-ens_list','--ens_list_file',help='File containing list of \
	ensemble pdb files to be used',type=str)
	parser.add_argument('-plot','--plot',help='Plot the distances?',type=bool)
	args=parser.parse_args()

	## Pre-processing inputs
	list_of_pdbs = args.ens_list_file.strip()

	try:
		res_range_HOX = map(int,args.hox_res_range.split('-'))
		res_range_flex = map(int,args.felx_res_range.split('-'))
	except ValueError as e:
		print 'Invalid input entered for -hox_range/-flex_range!'
		raise e

	f = open(args.ens_list_file,'r')
	print 'Reading list file: {}'.format(list_of_pdbs)
	list_lines = f.readlines()
	f.close()
	ens_pdb_names = [l.strip() for l in list_lines if l.strip().endswith('.pdb')]
	print 'Found {} pdbs in list file'.format(len(ens_pdb_names))
	
	distances_matrix = main(ens_pdb_names, res_range_HOX, res_range_flex)

	if args.plot:
		try:
			import matplotlib.pyplot as plt
			plot_results(distances_matrix)
		except ImportError:
			print 'Couldn"t import matplotlib for plotting!'

		f = open('raw_distances.json','w')
		f.write(json.dumps(distances_matrix))
		f.close()

		print 'Dumped the raw distances matrix to raw_distances.json'
