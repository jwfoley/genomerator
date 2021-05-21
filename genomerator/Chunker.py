import collections

class Chunker:
	'''
	given an iterable, accumulate its elements into chunks
	this class is itself an iterable, yielding chunks
	you define the breakpoint between chunks when instantiating this
	e.g. you can set a maximum length and the input will be broken into chunks of that length
	or you can say that the chunk ends if the next new feature is a certain distance from the last one in the chunk
	so it is not guaranteed that all yielded chunks will be the same length
	default chunk type is collections.deque since the only thing we know it will do for sure is append, but you can replace it with list or whatever else is more useful as long as it supports 'append'
	'''
	
	def __init__ (self,
		source, # iterable
		break_test, # function of current chunk and new feature that returns True if it is time to yield the current chunk and begin a new one
		chunk_type = collections.deque # type to use for chunks
	):
		self.source = source
		self.break_test = break_test
		self.chunk_type = chunk_type
		self.chunk = self.chunk_type()
		self._generator = self._yield_chunks()
	
	def _yield_chunks (self):
		for feature in self.source:
			if len(self.chunk) > 0 and self.break_test(self.chunk, feature):
				yield self.chunk
				self.chunk = self.chunk_type()
			self.chunk.append(feature)
		if len(self.chunk) > 0: yield self.chunk # could be empty if there was never a valid input
		self.chunk = None
	
	def __iter__ (self):
		return self
	
	def __next__ (self):
		return next(self._generator)

