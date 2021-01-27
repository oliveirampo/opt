import argparse

import action


def parse():
	parser = argparse.ArgumentParser()
	parser.add_argument("--type", dest="runType", required=True, type=str, choices=['GEN', 'ANA','OPT','PLOT'], help="Type of task to be done")
	parser.add_argument("--it", dest="it", required=True, type=int, help="Iteration")

	args = parser.parse_args()
	runType = args.runType
	it = args.it

	classes = {"GEN": action.Gen, "ANA": action.Ana, "OPT": action.Opt}
	act = classes[runType](it)
	return act

