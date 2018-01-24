import pdb_objects as pdb

### This is a very specific code. Better not to use it for general purposes.

def main(ens_pdb_names, res_range_flex):
	'''
	Takes care of pdb_parsing, center of mass and distance calculations
	'''

	all_calculations = []

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
		
		## Place holder for CA position vectors of all residues of flexible region
		flex_CA_vectors = []

		## Place holder for vectors joining one CA to next CA
		flex_bond_vectors = []

		## Loop through all the residues of chosen chain
		for residue in _chain:
			if res_range_flex[0]<=residue.serial<=res_range_flex[1]:
				flex_CA_vectors.append(residue.CA.position)
				# print 'Residue',residue.serial, 'is in flexible region.'

		## Make bond vectors by taking consecutive difference of position vectors.
		for i in range(len(flex_CA_vectors)-2):
			vec = flex_CA_vectors[i+1]-flex_CA_vectors[i]
			mod = abs(vec)
			flex_bond_vectors.append(vec/mod)

		calculations = {}

		no_of_CAs = len(flex_bond_vectors)-1

		for i in range(1,no_of_CAs+1):
			vleft, vright = flex_bond_vectors[i-1], flex_bond_vectors[i]
			s_range = range(1-i,no_of_CAs-i+1)
			temp = {}
			for s in s_range:
				if s<0:
					if abs(s)%2==1:
						temp[s] = vright.dot(flex_bond_vectors[i-1+s])
					else:
						temp[s] = vleft.dot(flex_bond_vectors[i-1+s])
				elif s==0:
					temp[s] = vleft.dot(vleft)
				elif s>0:
					if s%2==1:
						temp[s] = vleft.dot(flex_bond_vectors[i+s])
					else:
						temp[s] = vright.dot(flex_bond_vectors[i+s])

			calculations[i] = temp

		all_calculations.append(calculations)

	return all_calculations

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
	parser = argparse.ArgumentParser(description='CORRELATION FUNCTION CALCULATIONS FOR FLEXIBLE REGION')
	parser.add_argument("-flex_range",'--felx_res_range',help="Range of residues in pdb for flexible region")
	parser.add_argument('-ens_list','--ens_list_file',help='File containing list of \
	ensemble pdb files to be used',type=str)
	# parser.add_argument('-plot','--plot',help='Plot the distances?',type=bool)
	args=parser.parse_args()

	## Pre-processing inputs
	list_of_pdbs = args.ens_list_file.strip()

	try:
		res_range_flex = map(int,args.felx_res_range.split('-'))
	except ValueError as e:
		print 'Invalid input entered for -flex_range!'
		raise e

	f = open(args.ens_list_file,'r')
	print 'Reading list file: {}'.format(list_of_pdbs)
	list_lines = f.readlines()
	f.close()
	ens_pdb_names = [l.strip() for l in list_lines if l.strip().endswith('.pdb')]
	print 'Found {} pdbs in list file'.format(len(ens_pdb_names))
	
	all_calculations = main(ens_pdb_names, res_range_flex)

	f = open('correlation_test.json','w')

	avg_calculations = {}

	# for calc in all_calculations:
	# 	print calc

	for res in all_calculations[0]:
		temp = {}
		for s in all_calculations[0][res]:
			temp2 = 0
			for calc in all_calculations:
				temp2+=calc[res][s]
			temp[s] = temp2/len(all_calculations)
		avg_calculations[res] = temp

	f = open('correlation_beta_test.json','w')
	f.write(json.dumps(avg_calculations,indent=4,sort_keys=True))
	f.close()