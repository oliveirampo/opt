

class Atom:
	def __init__(self, idx, nam, iac, charge):
		self._idx = int(idx)
		self._nam = nam
		self._iac = int(iac)
		self._charge = charge
		self._ignore = False

	def __str__(self):
		s = '{:3} {:3} {} {} {}'.format(self._idx, self._iac, self._nam, self._bnd_nrm, self._charge.cur)
		return s

	@property
	def idx(self):
		return self._idx

	@property
	def nam(self):
		return self._nam

	@property
	def iac(self):
		return self._iac

	@property
	def charge(self):
		return self._charge

	@property
	def ignore(self):
		return self._ignore

	@ignore.setter
	def ignore(self, n):
		self._ignore = n
