import collections

class Buffer:
	'''
	given an iterable accumulate them into lists
	you define the breakpoint between lists when instantiating this
	e.g. you can set a size threshold to break input into roughly evenly sized chunks
	but you may also need to add a condition so you break in a safe place in the input
	so it is not guaranteed to return the same size chunks
	'''
	
	def __init__ (self,
		source, # iterable
		break_test # function to apply to self that determines when to yield a list of features
	):
		self.source = source
		self.break_test = break_test
		self.buffer = collections.deque()
		self._generator = self._yield_lists()
	
	def _yield_lists (self):
		for feature in self.source:
			self.buffer.append(feature)
			if self.break_test(self):
				yield self.buffer
				self.buffer = collections.deque()
		yield self.buffer
	
	def __iter__ (self):
		return self
	
	def __next__ (self):
		return next(self._generator)

