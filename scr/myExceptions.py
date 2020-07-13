

class MissingKeyWordError(KeyError):
	def __init__(self, key, src):
		self.key = key
		self.src = key
		s = "Missing key ({}) in {}".format(key, src)
		super().__init__(s)


class NoSuchFile(FileNotFoundError):
	def __init__(self, fileName):
		s = '\n\tNo such file: {}'.format(fileName)
		super().__init__(s)



