class Parameter():
	def __init__(self, typ, rng, ori, prev, cur):
		self._typ  = typ
		self._rng  = float(rng)
		self._ori  = float(ori)
		self._prev = float(prev)
		self._cur  = float(cur)


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



