class Matrix:
	def __init__(self):
		self._listIAC = {}

	@property
	def listIAC(self):
		return self._listIAC

	def add(self, iac1, iac2):
		iac1 = int(iac1)
		iac2 = int(iac2)

		if not iac1 in self._listIAC:
			self._listIAC[iac1] = [iac2]
		else:
			self._listIAC[iac1].append(iac2)


	def __str__(self):
		s = ''
		for iac in self._listIAC:
			s = s + '{} [{}]\n'.format(iac, ' '.join(map(str, self._listIAC[iac])))
		return s

	

