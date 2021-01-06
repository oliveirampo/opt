

class Atom:
	def __init__(self, idx, nam, iac, charge):
		self._idx = int(idx)
		self._nam = nam
		self._iac = int(iac)
		self._charge = charge
		self._bnd_nrm = {}
		self._bnd_nei = {}
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

	def addBndNrm(self, idx, dist):
		idx = int(idx)
		dist = float(dist)
		self._bnd_nrm.update({idx: dist})

	def addBndNei(self, idx, dist):
		idx = int(idx)
		dist = float(dist)
		self._bnd_nei.update({idx: dist})

	def isNrmNB(self, idx):
		idx = int(idx)
		if idx in self._bnd_nrm:
			return True
		return False

	def isNeiNB(self, idx):
		idx = int(idx)
		if idx in self._bnd_nei:
			return True
		return False

	def getNrmBndDist(self, idx):
		idx = int(idx)
		return self._bnd_nrm[idx]

	def getNeiBndDist(self, idx):
		idx = int(idx)
		return self._bnd_nei[idx]
