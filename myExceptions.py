class IOError(Exception):
	"""Base class for exceptions in this module."""
	pass


class MissingKeyWordError(IOError):
	def __init__(self, key, src):
		self.key = key
		self.src = src


