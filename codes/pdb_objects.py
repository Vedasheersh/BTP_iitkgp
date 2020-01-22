import numpy as np
from my_dict import OrderedDict

class PDBTree(object):
	"""
	Description:-
		Data structure to store the hierarchy of a .pdb file
	"""
	def __init__(self):
		super(PDBTree, self).__init__()
		self.models = OrderedDict()
		self.children = self.models
	
	def add(self,model_num):
		self.models[model_num] = Model(self,model_num)
		
	def model(self,model_num):
		'''
		Returns Model object.
		---------------------------------------------------------
		Args:
			`model_num` : `int` (model number)
		Usage:
			Eg: pdb_tree.model(1) or pdb_tree.model1
		'''
		try:
			return self.models[model_num]
		except KeyError as e:
			try:
				return self.models[str(model_num)]
			except KeyError as e:
				return None

	def __getattr__(self,arg):
		models_now = self.models.keys()
		if 'model' in arg and len(arg)>5:
			try:
				num = int(arg[5:])
				return self.model(num)
			except TypeError as e:
				return None
		if 'parent' in arg:
			return self

	def __str__(self):
		return 'PDBTree object: I have {} different model(s) as children.'.format(len(self.models))

	def __iter__(self):
		"""
		Utility:-
			Can iterate through self's children
		"""
		return self.children.values().__iter__()

## First level child object of PDBTree
class Model(object):
	"""
	Description:-
		child of PDBTree and parent to Chain(s)
	"""
	def __init__(self, pdbtree, model_num):
		"""
		Note:-
			Model can never be instantiated independently. It has to have
			a PDBTree instance as the parent.
		"""
		self.model_num = model_num
		self.chains = OrderedDict()
		if isinstance(pdbtree, PDBTree):
			self.parent = pdbtree
		self.children = self.chains

	def add(self, chain_name):
		self.chains[chain_name] = Chain(self, chain_name)
		return self.chains[chain_name]

	def get_dihedrals(self,chain_names=[]):
		dihedrals = {}
		if not chain_names:
			chain_names = self.chains.keys()
		for chain_name in chain_names:
			chain = self.chain(chain_name)
			dihedrals[chain_name] = chain.get_dihedrals()
		return dihedrals

	def chain(self,chain_name):
		'''
		Returns Chain object.
		---------------------------------------------------------
		Args:
			`chain_name` : `str` (chain identifier)
		Usage:
			TODO
		'''
		try:
			return self.chains[chain_name]
		except KeyError as e:
			return None

	def __getattr__(self,arg):
		models_now = self.chains.keys()
		if 'chain' in arg and len(arg)>5:
			try:
				chain_name = arg[5:]
				return self.chain(chain_name)
			except TypeError:
				return None
		
	def __iter__(self):
		return self.chains.values().__iter__()

	def __repr__(self):
		return 'Model object: Model {}'.format(self.model_num)

	def __str__(self):
		return 'Model object: I have {} chains as children'.format(len(self.children))
	
	def translate(self, xt, yt, zt, chains='ALL'):
        	if chains=='ALL':
            		chains = self.chains.keys()
        	for c in chains:
            		chain = self.chains[c]
            		chain.translate(xt,yt,zt)


class Chain(object):
	"""
	Description:-
		child of Model and parent to Residue(s)
	"""
	def __init__(self, model, chain_name):
		"""
		Note:-
			Chain can never be instantiated independently. It has to have
			a Model instance as the parent.
		"""
		self.chain_name = chain_name
		self.residues = OrderedDict()
		self.children = self.residues
		if isinstance(model, Model):
			self.parent = model		
			
	def add(self, res_num, res_name):
		self.residues[res_num] = Residue(self, res_num, res_name)
		return self.residues[res_num]

	def __iter__(self):
		return self.residues.values().__iter__()

	def residue(self, res_num):
		'''
		Returns Chain object.
		---------------------------------------------------------
		Args:
			`res_num` : `int` (residue number)
		Usage:
			TODO
		'''
		try:
			return self.residues[res_num]
		except KeyError as e:
			return None

	def get_dihedrals(self):
		dihedrals = []
		residues = self.residues.values()
		
		for res in residues:
			dihedrals.append([(res.res_num,res.res_name),res.get_dihedral()])

		return dihedrals

	def __getattr__(self,arg):
		models_now = self.residues.keys()
		if 'residue' in arg and len(arg)>7:
			try:
				res_num = int(arg[7:])
				return self.residue(res_num)
			except TypeError:
				#TODO
				return None
		elif 'res' in arg and len(arg)>3:
			try:
				res_num = int(arg[3:])
				return self.residue(res_num)
			except TypeError:
				#TODO
				return None

	def __str__(self):
		t = self.children.values()
		return 'Chain object: I have {} residues as children'.format(len(self.children))

	def __repr__(self):
		return 'Chain object: Chain {}'.format(self.chain_name)
	
	def translate(self, xt, yt, zt, residues='ALL'):
		if residues=='ALL':
			residues = self.residues.keys()
		for r in residues:
			residue = self.residues[r]
			residue.translate(xt,yt,zt)
		
class Residue(object):
	"""
	Description:-
		child of Chain and parent to Atom(s)
	"""
	def __init__(self, chain, res_num, res_name):
		"""
		Note:-
			Residue can never be instantiated independently. It has to have
			a Model instance as the parent.
		"""
		super(Residue, self).__init__()
		self.res_name = res_name
		self.res_num = res_num
		self.atoms = OrderedDict()
		self.children = self.atoms
		if isinstance(chain,Chain):
			self.parent = chain
			self.logger = self.parent.logger

	def add(self, serial, name):
		atom_name = name.strip()
		self.atoms[atom_name] = Atom(self,serial,name)
		return self.atoms[atom_name]

	def get_dihedral(self):
		try:
			prev_res = self.parent.residues[self.res_num-1]
		except KeyError:
			prev_res = None
		try:
			next_res = self.parent.residues[self.res_num+1]
		except KeyError:
			next_res = None

		N2vec = self.atoms['N'].position
		CA2vec = self.atoms['CA'].position
		C2vec = self.atoms['C'].position

		phi = 180
		psi = 180

		if prev_res:
			C1vec = prev_res.atoms['C'].position
			phi = C1vec.dihedral(N2vec,CA2vec,C2vec)
		
		if next_res:
			N3vec = next_res.atoms['N'].position
			psi = N2vec.dihedral(CA2vec,C2vec,N3vec)

		# print 'Phi:{}, Psi:{}'.format(phi,psi)
		return phi,psi

	def __iter__(self):
		return self.atoms.values().__iter__()

	def __str__(self):
		return 'Residue object: I have {} atoms as children'.format(len(self.children))

	def __repr__(self):
		return 'Residue object: Residue {} {}'.format(self.res_num, self.res_name)

	def __getattr__(self,arg):
		atom_names = ['N','CA','C','O']
		if arg in atom_names:
			return self.atoms[arg]
		elif arg=='name':
			return self.res_name
		elif arg=='serial':
			return self.res_num
	def translate(self, xt, yt, zt, atoms='ALL'):
		if atoms=='ALL':
			atoms = self.atoms.keys()
		for a in atoms:
			atom = self.get_atom(str(a))
			atom.translate(xt,yt,zt)

class Atom(object):
	"""
	Description:
		child of Residue and parent of None.
	"""
	def __init__(self, residue, serial, name):
		super(Atom, self).__init__()
		if isinstance(residue,Residue):
			self.parent= residue
		self.serial = serial
		self.name = name

	def add(self, coords, extras):
		if isinstance(coords, Vector):
			self.position = coords
		self.extras = extras
		return self.position

	def line(self):
		header = '{:<6}'.format('ATOM')
		serial = '{:>5}'.format(str(self.serial))+' '

		name = self.name+' '
		if len(name)==1 or len(name)==3:
			name = ' '+name
		name = '{:<4}'.format(str(name))
		
		res = self.parent.res_name+' '
		chain = self.parent.parent.chain_name

		res_num = '{:>4}'.format(str(self.parent.res_num))

		space = ' '*4

		coordinates = self.position.format_for_pdb()
		
		return header + serial + name + res + chain + res_num + space + coordinates + self.extras[:-1]

	def __repr__(self):
		return 'Atom object: {} {} in Residue {} {}'.format(self.serial,self.name,self.parent.res_num,self.parent.res_name)

	def __getattr__(self,arg):
		if arg=='name':
			return self.name
		if arg=='serial':
			return self.serial
	
	def translate( self, xt, yt, zt):
		tran_pos = Vector(xt,yt,zt)
		cur_pos = self.position
		self.position = cur_pos.translate(tran_pos)

class Vector(object):
	"""
	Description:
		Data structure to do vector algebra on coordinates (position) of atoms.
	"""
	def __init__(self, x=None,y=None,z=None):
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
		self.coords = (x,y,z)

	def __abs__(self):
		return np.sqrt(self.x**2+self.y**2+self.z**2)

	def __str__(self):
		return 'Vector:<{}i {}j {}k>'.format(self.x,self.y,self.z)

	def __add__(self,other):
		return Vector(self.x+other.x,self.y+other.y,other.z+self.z)

	def __sub__(self,other):
		return Vector(self.x-other.x,self.y-other.y,self.z-other.z)

	def __mul__(self,other):
		if isinstance(other,int) or isinstance(other,float):
			x = self.x*other
			y = self.y*other
			z = self.z*other
		return Vector(x,y,z)

	def __div__(self, other):
		if isinstance(other,int) or isinstance(other,float):
			x = self.x/other
			y = self.y/other
			z = self.z/other
		return Vector(x,y,z)

	def alter(self,other):
		self.__init__(*other.coords)

	def dot(self,other):
		return np.dot(np.array(self.coords),np.array(other.coords))

	def cross(self,other):
		return Vector(*(np.cross(np.array(self.coords),np.array(other.coords))))

	def angle(self,other):
		return 180*np.arccos(self.dot(other)/abs(self)/abs(other))/np.pi

	def distance(self,other):
		return abs(self-other)

	def dihedral(a,b,c,d):
		ab = b-a
		bc = c-b
		cd = d-c
		norm1 = ab.cross(bc)
		norm2 = bc.cross(cd)
		angle = norm1.angle(norm2)
		foo=bc.dot(norm1.cross(norm2))
		# print foo
		if foo>0:
			return angle
		elif foo<0:
			return -1*angle

	def translate(self,translate_by):
		'''
		Returns translated coordiates of self
		'''
		try:
			return self-translate_by
		except:
			print 'translate_by has to be a Coord object as well!'

	def rotate(self,rotate_by):
		'''
		Args:-
			rotate_by: 3x3 np array or np matrix containing the rotational operator.
		'''
		a = np.matrix(self.coords)
		r = np.matrix(rotate_by)

		b = np.dot(a,r)
		return Vector(*b.tolist()[0])

	def format_for_pdb(self):
		rounded = (round(self.x,3),round(self.y,3),round(self.z,3))
		str_rounded = map(str,rounded)
		f_str = '{:>8}'
		return f_str.format(str_rounded[0])+f_str.format(str_rounded[1])+f_str.format(str_rounded[2])

class PDBParser(object):
	"""
	Parse a .pdb file and make its PDBTree.
	"""
	def __init__(self,pdb_file=None):
		'''
		Create a parser object.
		'''
		super(PDBParser, self).__init__()
		
		if pdb_file:
			parse(pdb_file)

	def parse(self,pdb_file):
		'''
		Parse a .pdb file.
		---------------------------------------------
			Args:
				pdb_file = `str` - pdb file name
		'''
		def _line_reader():
			'''
			Generates one line at a time from the file.
			'''
			for line in file_handle:
				yield line

		def _parse_line(line):
			'''
			Parse a single line from the file.
			'''
			try:
				serial, atom_name, res_name, chain_id, res_num, coords, extras = \
				line[6:11], line[12:16], line[17:20], line[21], line[22:26], (line[30:38],line[38:46],line[46:54]), line[54:]
				return serial, atom_name, res_name, chain_id, res_num, coords, extras
			except IndexError:
				#TODO
				pass

		## Open the pdb file
		try:
			file_handle = open(pdb_file)
			print 'Opened pdb file: {}'.format(pdb_file)
		except IOError as e:
			print 'Pdb file: {} not found!'.format(pdb_file)

		reader = _line_reader()

		line_num = 0
		file_end = False

		# Initialize PDBTree object with logger
		pdb = PDBTree()
		# Add the first model 
		pdb.add(0)
		cur_model = pdb.model(0)

		## To find first ATOM record and instantiate objects
		while True:
			try:
				line = next(reader)
				line_num += 1
				if line.startswith('ATOM'):
					_serial, _atom_name, _res_name, _chain_id, _res_num, _coords,_extras = _parse_line(line)
					_chain = cur_model.add(_chain_id)
					_res = _chain.add(int(_res_num), _res_name)
					_atom = _res.add(_serial, _atom_name)
					_atom.add(Vector(*_coords),_extras)
					print 'First ATOM record found in line {}'.format(line_num)
					## Only one ATOM record is read!
					break

			except StopIteration:
				file_end = line_num
				print 'No ATOM records found in the pdb!'
				break
		
		while not file_end:
			## Next line from the file
			try:
				line = next(reader)
				line_num += 1

				if line.startswith('ATOM'): 

					serial, atom_name, res_name, chain_id, res_num, coords, extras = _parse_line(line)

					if chain_id==_chain_id and res_num==_res_num:
						## No chain/res change
						atom = _res.add(serial,atom_name)
						atom.add(Vector(*coords),extras)

					elif chain_id==_chain_id and res_num!=_res_num:
						## No chain change but different residue
						_res = _chain.add(int(res_num), res_name)
						atom = _res.add(serial,atom_name)
						atom.add(Vector(*coords),extras)

					elif chain_id!=_chain_id:
						## Different chain so, different residue also
						_chain = cur_model.add(chain_id)
						_res = _chain.add(int(res_num),res_name)
						atom = _res.add(serial,atom_name)
						atom.add(Vector(*coords),extras)
						
					_chain_id = chain_id
					_res_num = res_num
					
			except StopIteration:
				file_end = line_num
				break
		
		return pdb
