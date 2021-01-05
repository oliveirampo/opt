class Parameter:
	def __init__(self, iac, typ, rng, ori, prev, cur):
		self._iac = iac
		self._typ = typ
		self._rng = float(rng)
		self._ori = float(ori)
		self._prev = float(prev)
		self._cur = float(cur)
		self._min = float(ori)
		self._max = float(ori)
		self._symmetry = ''

	@property
	def iac(self):
		return self._iac

	@property
	def typ(self):
		return self._typ

	@property
	def rng(self):
		return self._rng

	@property
	def ori(self):
		return self._ori

	@property
	def prev(self):
		return self._prev

	@property
	def cur(self):
		return self._cur

	@property
	def min(self):
		return self._min

	@property
	def max(self):
		return self._max

	@property
	def symmetry(self):
		return self._symmetry

	@cur.setter
	def cur(self, val):
		self._cur = val

	@min.setter
	def min(self, val):
		self._min = val

	@max.setter
	def max(self, val):
		self._max = val

	@symmetry.setter
	def symmetry(self, val):
		self._symmetry = val

	def hasSymmetricIAC(self):
		return self._symmetry != ''
