"""Module with exceptions.

Classes:
	MissingKeyWordError
	NoSuchFile
	ClassNotImplemented
"""


class MissingKeyWordError(KeyError):
	"""Mapping key not found."""

	def __init__(self, key, src):
		self.key = key
		self.src = key
		s = "Missing key ({}) in {}".format(key, src)
		super().__init__(s)


class NoSuchFile(FileNotFoundError):
	"""File not found."""
	def __init__(self, fileName):
		s = '\n\tNo such file: {}'.format(fileName)
		super().__init__(s)


class ClassNotImplemented(NotImplementedError):
	"""Method or function hasn't been implemented yet."""
	def __init__(self, className, functionName):
		s = '\n\tClass not implemented ({}) in module {}'.format(className, functionName)
		super().__init__(s)
