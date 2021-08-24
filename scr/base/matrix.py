"""Matrix that defines if C12(II) parameters is used (in case of hydrogen-bonding interactions).

Class:
	Matrix
"""


class Matrix:
	"""Matrix that defines if C12(II) parameters is used
	(in case of hydrogen-bonding interactions).

	Attributes:
		_listIAC: (dict) Dictionary where
		key = atom type code
		value = list of atoms with which the C12(II) parameter are used.

	Methods:
		add():

	"""

	def __init__(self):
		"""Constructs all the necessary attributes for the matrix."""

		self._listIAC = {}

	@property
	def listIAC(self):
		return self._listIAC

	def add(self, iac1, iac2):
		"""Adds iac2 to list of atom types of iac1.

		:param iac1: (int)
		:param iac2: (int)
		:return:
		"""

		iac1 = int(iac1)
		iac2 = int(iac2)

		if iac1 not in self._listIAC:
			self._listIAC[iac1] = [iac2]
		else:
			self._listIAC[iac1].append(iac2)

	def __str__(self):
		s = ''
		for iac in self._listIAC:
			s = s + '{} [{}]\n'.format(iac, ' '.join(map(str, self._listIAC[iac])))
		return s

	def __contains__(self, item):
		"""Allows the following type of comparison:
			if iac in matrx: ...

		:param item: (int) Atom type index.
		:return: (boolean)
		"""
		return item in self._listIAC

	def __getitem__(self, item):
		"""Returns list associated with a given atom type index.

		:param item: (int) Atom type index.
		:return: (list) List of atom types.
		"""
		return self._listIAC[item]
