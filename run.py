import sys

from configuration import Conf
import myExceptions
import inpParser


def main():
	try:
		action = inpParser.parse()
		conf = Conf('00_inp/Conf.dat', action.it)

		molecules, atomTypes = action.read_inp_files(conf)
		action.run(conf, molecules, atomTypes)

	except (FileNotFoundError) as err:
		print(err)
		sys.exit('\n')

	except myExceptions.MissingKeyWordError as err:
		print("\n\tNo keyword ({}) in {}".format(err.key, err.src))
		sys.exit('\n')


if __name__ == "__main__":
	main()