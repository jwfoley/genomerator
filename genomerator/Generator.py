import collections
from .GenomeFeature import GenomeFeature

class RegionGenerator (object):
	'''
	given a list of reference lengths, generates GenomeFeatures representing regions of a specified length
	'''
	
	def __init__ (self,
		reference_lengths, 
		region_length =    1, # length of each region
		overlap  =         False, # return consecutive overlapping regions offset by 1 base
		include_partial =  True, # return shorter region at the end of the reference if not exactly divisible
		default_factory =  int # function to create the default .data of each region
	):
		self.reference_lengths = reference_lengths
		self.region_length = region_length
		self.overlap = overlap
		self.include_partial = include_partial
		self.default_factory = default_factory
		self._generator = self._yield_regions()
	
	def _yield_regions (self):
		for reference_id, reference_length in zip(range(len(self.reference_lengths)), self.reference_lengths):
			for left_pos in range(1, reference_length + 1, (1 if self.overlap else self.region_length)):
				if left_pos + self.region_length - 1 > reference_length: # new feature would extend past end of reference
					if self.include_partial:
						right_pos = reference_length
					else:
						break
				else:
					right_pos = left_pos + self.region_length - 1
				yield GenomeFeature(
					reference_id =  reference_id,
					left_pos =      left_pos,
					right_pos =     right_pos,
					data =          self.default_factory()
				)
	
	def __iter__ (self):
		return self
	
	def __next__ (self):
		return next(self._generator)


class OperationGenerator (object):
	'''
	for every feature in iterable "b", perform some operation on the corresponding feature(s) in iterable "a", then yield the altered features from iterable "a"
	e.g. "a" might be a list of genome regions and "b" might be a list of sequence reads, and then you count the number of reads intersecting each region
	and that is what the default configuration will do
	but you can customize the functions to test matches, e.g. maybe you only care if the "b" feature is entirely inside the "a", or only its first base is
	you can also customize the functions that decide when to stop checking more features, e.g. maybe you want to allow hits some distance past the end of a region
	this probably only makes sense if both iterables are sorted, but doesn't actually check
	compare bedtools
	'''
	
	def __init__ (self,
		a, # iterable of GenomeFeature instances to yield after modifying
		b, # iterable of GenomeFeature instances to use for modifying "a"
		operate =                lambda a,b: a.__setattr__('data', a.data + 1), # function to apply to a matching feature from "a" given a new feature from "b"
		match =                  lambda a,b: a.intersects(b), # function to check whether a feature from "a" matches a feature from "b"
		check_a_passed =         lambda a,b: a.left_of(b), # function to check whether a feature from "a" is done being updated (e.g. we've completely passed it), given the latest feature from "b"
		check_b_passed =         lambda a,b: a.right_of(b), # function to check whether a feature from "b" is done finding matches (e.g. we've completely passed it), given the latest feature from "a" 
		stop_at_first_match =    False, # for each feature from "b", only perform the operation on a single feature from "a", then stop looking for more matches (don't double-count)
	):
		self.a = a
		self.b = b
		self.operate = operate
		self.match = match
		self.check_a_passed = check_a_passed
		self.check_b_passed = check_b_passed
		self.stop_at_first_match = stop_at_first_match
		self._a_features = collections.deque()
		self._generator = self._yield_features()
	
	def _yield_features (self):
		for b_feature in self.b:
			# first, yield any "a" features that are now done
			while len(self._a_features) > 0 and self.check_a_passed(self._a_features[0], b_feature):
				yield self._a_features.popleft()
			
			# second, add all "a" features necessary to cover this "b" feature (plus one more)
			if len(self._a_features) == 0 or not self.check_b_passed(self._a_features[-1], b_feature):
				for a_feature in self.a:
					if self.check_a_passed(a_feature, b_feature): # don't bother with "a" features "b" has already passed
						yield a_feature
					else:
						self._a_features.append(a_feature)
						if self.check_b_passed(a_feature, b_feature): break
			
			# third, update the current "a" features according to this "b" feature
			for a_feature in self._a_features:
				if self.match(a_feature, b_feature): # found a match
					self.operate(a_feature, b_feature)
					if self.stop_at_first_match: break # stop looking
				elif self.check_b_passed(a_feature, b_feature): # we've already passed the "b" feature
					break # so stop looking
		
		# purge all remaining "a" features because there is no more "b"
		for a_feature in self._a_features: yield a_feature
	
	def __iter__ (self):
		return self
	
	def __next__ (self):
		return next(self._generator)

