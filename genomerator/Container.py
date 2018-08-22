import collections
from .GenomeFeature import GenomeFeature

class GenomePositionDeque (GenomeFeature):
	'''
	deque wrapper in which each element corresponds to a genome position
	if you add items to the ends or remove them from the ends, the coordinates change appropriately
	this can mean extending past the ends of the reference sequence!
	inserting, removing, etc. is not allowed because it's ambiguous which way to move the coordinates
	'''
	
	def __init__ (self, *arg, **kwarg):
		super().__init__(*arg, **kwarg)
		if self.data is not None:
			assert len(self.data) == len(self)
			self.data = collections.deque(self.data)
		else:
			self.data = collections.deque([None] * len(self))	

	# position lookup
	
	def __getitem__ (self, index):
		return self.data[index]
	
	def __setitem__ (self, index, value):
		self.data[index] = value
	
	def get_pos (self, position):
		'''
		look up a value by genome position rather than list position
		'''
		if not self.left_pos <= position <= self.right_pos: raise IndexError
		return self.data[position - self.left_pos]
	
	def set_pos (self, position, value):
		'''
		set a value by genome position rather than list position
		'''
		if not self.left_pos <= position <= self.right_pos: raise IndexError
		self.data[position - self.left_pos] = value
	
	def index (self, value, start = None, stop = None):
		if start is None: start = 0
		if stop is None: stop = len(self)
		return self.data.index(value, start, stop)
	
	def position (self, value, start = None, stop = None):
		'''
		like index but for a genome position
		note that for consistency with index, stop is the position *after* the end of the range
		'''
		if start is None: start = self.left_pos
		if stop is None: stop = self.right_pos + 1 # because it needs to be past the end
		if not (self.left_pos <= start <= self.right_pos and self.left_pos <= stop <= self.right_pos + 1): raise IndexError
		return self.data.index(value, start - self.left_pos, stop - self.left_pos) + self.left_pos
	
	def __iter__ (self):  
		'''
		if you iterate over this class, it produces GenomeFeatures each containing the data from one genome position
		iterate over self.data if you just want the deque contents
		'''
		return (GenomeFeature(
			reference_id =  self.reference_id,
			left_pos =      index + self.left_pos,
			right_pos =     index + self.left_pos,
			data =          value
		) for index, value in zip(range(len(self)), self.data))
		
	
	# adding position data
	
	def append (self, value):
		self.data.append(value)
		self.right_pos += 1
	
	def appendleft (self, value):
		self.data.appendleft(value)
		self.left_pos -= 1
	
	def extend (self, values):
		self.data.extend(values)
		self.right_pos += len(values)
	
	def extendleft (self, values):
		self.data.extendleft(values)
		self.left_pos -= len(values)
	
	
	# removing position data
	
	def pop (self):
		self.right_pos -= 1
		return self.data.pop()
	
	def popleft (self):
		self.left_pos += 1
		return self.data.popleft()
		
	
	# full-length operations
	
	def count (self, value):
		return self.data.count(value)
	
	def reverse (self):
		return self.data.reverse()
	
	def rotate (self, steps = 1):
		self.data.rotate(steps)
		self.left_pos += steps
		self.right_pos += steps
	
	def copy (self):
		return self.__class__(
			reference_id =  self.reference_id,
			left_pos =      self.left_pos,
			right_pos =     self.right_pos,
			data =          self.data.copy()
		)
	
	def intersection(self, other):
		'''
		return an instance of the overlap between this deque and another genome region, including only the data in the intersection, or None if no overlap
		use this to slice part of the original region too
		'''
		if other.reference_id != self.reference_id: return None
		if other.right_pos < self.left_pos or self.right_pos < other.left_pos: return None
		new_left = max(self.left_pos, other.left_pos)
		new_right = min(self.right_pos, other.right_pos)
		return self.__class__(
			reference_id =  self.reference_id,
			left_pos =      new_left,
			right_pos =     new_right,
			data =          collections.deque(self.data[i] for i in range(new_left - self.left_pos, new_right - self.left_pos)) # remember, deques can't be directly sliced
		)


class GenomeBuffer (GenomeFeature):
	'''
	a generator that maintains an internal buffer (deque) of coordinate-sorted GenomeFeatures, harvested from some other iterable, and pops them off when some condition is met (e.g. newest feature doesn't overlap the oldest, implying no future features will overlap it either)
	left and right pos track the range of *left* positions of all features in the buffer (some of them may extend outside it)
	'''
	
	__slots__ = '_feature_iterable'
	
	def __init__ (self, feature_iterable):
		self._feature_iterable = feature_iterable
		first_feature = next(self._feature_iterable)
		super().__init__(
			reference_id =  first_feature.reference_id,
			left_pos =      first_feature.left_pos,
			right_pos =     first_feature.right_pos,
			data =          collections.deque([first_feature])
		)
	
	def advance (self):
		try:
			next_feature = next(self._feature_iterable)
			if len(self.data) > 0 and not next_feature >= self.data[-1]: raise RuntimeError('features not sorted correctly')
		except StopIteration:
			next_feature = None

		output = []
		while len(self.data) > 0 and (next_feature is None or next_feature.right_of(self)):
			output.append(self.data.popleft())
		
		if next_feature is not None:
			self.data.append(next_feature)
			self.reference_id =  next_feature.reference_id
			self.left_pos =      self.data[0].left_pos
			self.right_pos =     max(self.right, next_feature.right).right_pos # do the comparison on GenomeFeatures so it takes reference_id into account (don't use a higher position from a lower reference_id)
			return output
		else:
			raise StopIteration

