class OrderedDict(dict):
	"""
	Description:
		A wrapper object around dict which stores additional info. of the order in which (key,value) pairs are added to it.
	"""
	def __init__(self):
		self.data = {}
		self.order = []

	def __setitem__(self, key, item): 
		self.data[key] = item
		self.order.append(key)

	def __delitem__(self, key): 
		del self.data[key]
		self.order.remove(key)

	def __iter__(self):
		return self.order.__iter__()

	def next(self):
		"""
		Note:-
			This 'self' is the return value of __iter__
		"""
		for key in self:
			return self.data[key]

	def __repr__(self): return repr(self.data)
	__hash__ = None # Avoid Py3k warning
	def __len__(self): return len(self.data)
	def __getitem__(self, key):
		if key in self.data:
			return self.data[key]
		if hasattr(self.__class__, "__missing__"):
			return self.__class__.__missing__(self, key)
		raise KeyError(key)
	def keys(self): return [key for key in self.order]
	def items(self): return [(key,self.data[key]) for key in self.order]
	def values(self): return [self.data[key] for key in self.order]
	def has_key(self, key): return key in self.order
	def __contains__(self, key): return key in self.order
