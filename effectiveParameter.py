class EffectiveParameter:
	def __init__(self, idx, typ, iac1, iac2, typAtm1, typAtm2, val):
		self._idx     = int(idx)
		self._typ     = typ
		self._iac1    = int(iac1)
		self._iac2    = int(iac2)
		self._typAtm1 = typAtm1
		self._typAtm2 = typAtm2
		self._oriVal  = val
		self._prevVal = val
		self._curVal  = val
		self._used    = False

	@property
	def idx(self):
		return self._idx
	@property
	def typ(self):
		return self._typ
	@property
	def iac1(self):
		return self._iac1
	@property
	def iac2(self):
		return self._iac2
	@property
	def typAtm1(self):
		return self._typAtm1
	@property
	def typAtm2(self):
		return self._typAtm2
	@property
	def oriVal(self):
		return self._oriVal
	@property
	def prevVal(self):
		return self._prevVal
	@property
	def curVal(self):
		return self._curVal
	@property
	def used(self):
		return self._used

# TODO - convert this file to atom and LJPair


class CHG(EffectiveParameter):
	def __init__(self, idx, typ, iac, typAtm, val):
		EffectiveParameter.__init__(self, idx, "CHG_ATM", 1, iac, "MOLEC", typAtm, val)


class LJPair(EffectiveParameter):
	def __init__(self, idx, typ, iac1, iac2, typAtm1, typAtm2, val):
		EffectiveParameter.__init__(self, idx, typ, iac1, iac2, typAtm1, typAtm2, val)

