"""Atom class.

Classes:
--------
	Atom
"""

class Atom:
	"""Atom class."""

	def __init__(self, idx, nam, iac, charge):
		"""Constructs all the necessary attributes for the given atom.

		:param idx: (int) Index.
		:param nam: (str) Name
		:param iac: (int) Atom type index.
		:param charge: (Charge) Charge of atom.
		"""

		self._idx = int(idx)
		self._nam = nam
		self._iac = int(iac)
		self._charge = charge
		self._ignore = False

	def __str__(self):
		s = '{:3} {:3} {} {}'.format(self._idx, self._iac, self._nam, self._charge.cur)
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
