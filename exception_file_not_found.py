import logging

def get_number():
	return int('foo')

try:
	x = get_number()
except Exception as ex:
	logging.excetion('Caught an error')
