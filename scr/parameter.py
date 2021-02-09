"""Module for Parameter object.

Classes:
	Parameter
"""

class Parameter:
	"""Defines a parameter object.

	Attributes:
		_iac: (int) Atom type index.
		_typ: (float) Type of parameter.
		_rng: (float) Range of parameter variation.
		_ori: (float) Original value of parameter.
		_prev: (float) Previous value of parameter.
		_cur: (float) Current value of parameter.
		_min: (float) Minimal value of parameter.
		_max: (float) Maximal value of parameter.
		_symmetry = (str) Symmetry (code) of parameter.

	Methods:
		hasSymmetricIAC()
	"""

	def __init__(self, iac, typ, rng, ori, prev, cur):
		"""Constructs all the necessary attributes for the parameter.

		:param iac: (int) Atom type index.
		:param typ: (float) Type of parameter.
		:param rng: (float) Range of parameter variation.
		:param ori: (float) Original value of parameter.
		:param prev: (float) Previous value of parameter.
		:param cur: (float) Current value of parameter.
		"""

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
		"""Checks if parameters has a symmetric parameter."""

		return self._symmetry != ''
