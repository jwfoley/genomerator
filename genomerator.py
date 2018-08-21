import collections

# look at scipy.integrate.ode
# make a wiggle reader
# make a FASTA reader: optionally set span (None for just whatever the line length), overlap (if span > 1, should features overlap or just be contiguous, and if overlap how much?), behavior at end of each reference (if span > 1, return a final truncated feature?)
# obviously this also means there should be a FASTA writer (just set line length?)
# this also suggests 2bit should be supported
# still looking for sequential-read bigBed/bigWig module


class GenomeFeature (object):
	'''
	container for a generic genome feature with reference ID number (not name), left and right positions (1-based), strand (as 'is_reverse'), and embedded data of any kind
	if no right_pos is specified, it represents a single genome position (and right_pos = left_pos)
	non-stranded features default to the forward strand (is_reverse = False)
	is it possible to get to self.data by dereferencing the instance like a C++ iterator? does this have to be a subclass of a pointer or something?
	'''
	
	__slots__ = 'reference_id', 'left_pos', 'right_pos', 'is_reverse', 'data'
	
	# constructors
	
	def __init__ (self,
		reference_id,        # index number of the reference (chromosome), 0-based
		left_pos,            # leftmost position, 1-based (!)
		right_pos = None,    # rightmost position, 1-based
		is_reverse = False,  # True if on reverse strand, False if forward strand or unstranded
		data = None          # optional data (can be whatever object you want)
	):
		assert isinstance(reference_id, int) and reference_id >= 0
		assert isinstance(left_pos, int)
		self.reference_id = reference_id
		self.left_pos = left_pos
		if right_pos is None:
			self.right_pos = left_pos
		else:
			assert isinstance(right_pos, int)
			assert(right_pos >= left_pos)
			self.right_pos = right_pos
		self.is_reverse = is_reverse
		self.data = data
	
	@classmethod
	def from_genomefeature (cls, other, data = None):
		'''
		initiate from another GenomeFeature but optionally change the data
		'''
		return cls(
			reference_id =  other.reference_id,
			left_pos =      other.left_pos,
			right_pos =     other.right_pos,
			is_reverse =    other.is_reverse,
			data =          data if data is not None else other.data
		)
	
	@classmethod
	def from_alignedsegment (cls, aligned_segment):
		'''
		convert a pysam.AlignedSegment into a GenomeFeature
		the original AlignedSegment is still stored as self.data
		'''
		return cls(
			reference_id =  aligned_segment.reference_id,
			left_pos =      aligned_segment.reference_start + 1,
			right_pos =     aligned_segment.reference_end,
			is_reverse =    aligned_segment.is_reverse,
			data =          aligned_segment
		)
	
	@classmethod
	def from_bed (cls, bed_line, references):
		'''
		make a GenomeFeature from a line of a BED file
		you must give a list of reference names, in order, so it can find the index
		returns the split but unparsed fields in self.data
		consider making data an ordered dictionary for easy un-parsing
		but then also consider allowing user to specific the names of extra nonstandard fields (this will imply the number of extra fields, and then when parsing the file we will know the number of standard fields after splitting and subtracting that)
		''' 
		fields = bed_line.rstrip().split()
		return cls(
			reference_id =  references.index(fields[0]),
			left_pos =      int(fields[1]) + 1,
			right_pos =     int(fields[2]),
			is_reverse =    len(fields) >= 6 and fields[5] == '-',
			data =          fields
		)
	
	
	# accessing the attributes
	
	@property
	def left (self):
		'''
		return a new instance containing only the leftmost position, for quick comparisons
		'''
		return GenomeFeature(reference_id = self.reference_id, left_pos = self.left_pos, is_reverse = self.is_reverse, data = self.data) # should this be self.__class__? the problem is that it doesn't work in inheritance
	
	@property
	def right (self):
		'''
		return a new instance containing only the rightmost position, for quick comparisons
		'''
		return GenomeFeature(reference_id = self.reference_id, left_pos = self.right_pos, is_reverse = self.is_reverse, data = self.data) # should this be self.__class__? the problem is that it doesn't work in inheritance
	
	@property
	def start (self):
		'''
		return a new instance containing only the start position (left if forward strand, right if reverse strand)
		'''
		return self.right if self.is_reverse else self.left
	
	@property
	def end (self):
		'''
		return a new instance containing only the end position (right if forward strand, left if reverse strand)
		'''
		return self.left if self.is_reverse else self.right
	
	@property
	def start_pos (self):
		'''
		return the start position as a single coordinate
		'''
		return self.right_pos if self.is_reverse else self.left_pos
	
	@property
	def end_pos (self):
		'''
		return the end position as a single coordinate
		'''
		return self.left_pos if self.is_reverse else self.right_pos	
	
	# consider adding a __getitem__ by index i that returns a new instance containing the ith position in the range
	
	def __len__ (self):
		return self.right_pos - self.left_pos + 1
	
	
	# modifying the coordinates (and returning a new modified instance)
	
	def shift_left (self, distance):
		'''
		return a new instance with the left coordinate shifted the specified distance
		'''
		return GenomeFeature(reference_id = self.reference_id, left_pos = self.left_pos + distance, right_pos = self.right_pos, is_reverse = self.is_reverse, data = self.data)
	
	def shift_right (self, distance):
		'''
		return a new instance with the right coordinate shifted the specified distance
		'''
		return GenomeFeature(reference_id = self.reference_id, left_pos = self.left_pos, right_pos = self.right_pos + distance, is_reverse = self.is_reverse, data = self.data)
	
	# consider something for shifting forward/backward depending on feature strand
	
	def __add__ (self, distance):
		'''
		return a new instance with both coordinates shifted the specified distance
		'''
		return GenomeFeature(reference_id = self.reference_id, left_pos = self.left_pos + distance, right_pos = self.right_pos + distance, is_reverse = self.is_reverse, data = self.data)

	def __sub__ (self, distance):
		'''
		opposite of __add__
		'''
		return self.__add__(-distance)
	
	def switch_strand (self):
		'''
		return a new instance with the opposite strandedness
		'''
		return GenomeFeature(reference_id = self.reference_id, left_pos = self.left_pos, right_pos = self.right_pos, is_reverse = not self.is_reverse, data = self.data)
	
	
	# comparisons
	
	def __eq__ (self, other): # note right_pos is not checked (for sort order)
		return self.reference_id == other.reference_id and self.left_pos == other.left_pos
	
	def __lt__ (self, other): # note right_pos is not checked (for sort order)
		return (
			self.reference_id < other.reference_id or (
				self.reference_id == other.reference_id and
				self.left_pos < other.left_pos
			)
		)
	
	def __le__ (self, other):
		return self == other or self < other
	
	def left_of(self, other):
		'''
		this feature is completely to the left of the other with no overlap
		(or on an earlier reference altogether)
		'''
		return self.right < other.left
	
	def right_of (self, other):
		'''
		this feature is completely to the right of the other with no overlap
		(or on a later reference altogether)
		'''
		return self.left > other.right
	
	def same_as (self, other):
		'''
		this feature has the same start and end positions and strand as the other
		'''
		return self.left == other.left and self.right == other.right and self.is_reverse == other.is_reverse
	
	def intersection (self, other):
		'''
		return a GenomeFeature of the overlap between two regions (preserves this one's data), or None if no overlap
		check this for bugs!
		'''
		if other.reference_id != self.reference_id: return None
		if other.right_pos < self.left_pos or self.right_pos < other.left_pos: return None
		return GenomeFeature( # should this be self.__class__? the problem is that it doesn't work in inheritance
			reference_id =  self.reference_id,
			left_pos =      max(self.left_pos, other.left_pos),
			right_pos =     min(self.right_pos, other.right_pos),
			data =          self.data
		)
	
	def union (self, other):
		'''
		return a GenomeFeature of the combined region spanned by two GenomeFeatures
		returns None if they are on different references or there is a gap between them
		'''
		if (
			other.reference_id != self.reference_id or
			other.right_of(self.shift_right(1)) or # shift to account for two features that don't overlap but have no gap between them
			other.left_of(self.shift_left(-1))
		):
			return None
		else:
			return GenomeFeature (
				reference_id =  self.reference_id,
				left_pos =      min(self.left_pos, other.left_pos),
				right_pos =     max(self.right_pos, other.right_pos),
				data =          self.data
			)
	
	
	# progress (proportion of genome traversed)
	# make a class of progress tracker that stores the list
	
	def progress (self, reference_lengths):
		total_length = sum(reference_lengths)
		assert total_length > 0
		traversed_length = sum(reference_lengths[:self.reference_id]) + self.left_pos
		return traversed_length / total_length
	
	def progress_percent (self, reference_lengths):
		return round(100 * self.progress(reference_lengths))
	
	
	# exporting
	# use a class to store the names
	
	def bed (self,
		reference_names,
		name =          None,
		score =         None,
		thick_range =   None, # both 1-indexed
		rgb =           None,
		block_count =   None,
		block_sizes =   None,
		block_starts =  None, # 0-indexed relative to the leftmost position, as usual
		use_strand =    True
	):
		'''
		format the feature as a UCSC BED line
		requires a list of reference names since the feature only knows its index
		only returns reference name, coordinates, and strand by default, but you can specify the other fields or disable the strand
		'''
	
		assert thick_range is None or len(thick_range) == 2
		assert block_count is None or len(block_sizes) == len(block_starts) == block_count
		assert rgb is None or len(rgb) == 3
		
		# first make a complete entry
		fields = [
			reference_names[self.reference_id],                                    # chrom
			str(self.left_pos - 1),                                                # chromStart
			str(self.right_pos),                                                   # chromEnd
			str(name) if name is not None else '.',                                # name
			str(score) if score is not None else '.',                              # score
			('-' if self.is_reverse else '+') if use_strand else '.',              # strand
			str(thick_range[0] + 1) if thick_range is not None else '.',           # thickStart
			str(thick_range[1]) if thick_range is not None else '.',               # thickEnd
			','.join(map(str, rgb)) if rgb is not None else '.',                   # itemRgb
			str(block_count) if block_count is not None else '.',                  # blockCount
			','.join(map(str, block_sizes)) if block_sizes is not None else '.',   # blockSizes
			','.join(map(str, block_starts)) if block_starts is not None else '.'  # blockStarts
		]
		# now trim off the empty fields
		final_length = len(fields)
		for field in reversed(fields):
			if field == '.':
				final_length -= 1
			else:
				break
		
		return '\t'.join(fields[:final_length])
	
	def bedgraph (self, reference_names):
		'''
		format the feature as a UCSC BedGraph line (https://genome.ucsc.edu/goldenPath/help/bedgraph.html)
		remember to create a track header line!
		outputs whatever is in self.data as the dataValue; be sure this is something numerical because this method doesn't check
		'''
		
		return '\t'.join(reference_names[self.reference_id], str(self.left_pos - 1), str(self.data))
	
	def __repr__ (self):
		return ('%s(reference_id = %i, left_pos = %i, right_pos = %i, is_reverse = %s, data = %s)' % (self.__class__.__name__, self.reference_id, self.left_pos, self.right_pos, self.is_reverse, self.data))


class FeatureStream (object):
	'''
	given an iterable of GenomeFeature instances, yield them back
	but also verify correct sort order, if desired, and keep a count
	optionally replace the data with something user-specified instea
	optionally filter the data to only includ
	what if this could be a subclass of GenomeFeature so its position etc. can be easily compared? risky but interesting
	consider writing this so it just knows which type of data is coming in instead of subclasses
	'''
	
	__slots__ = 'source', 'default_data', 'filter', 'verify_sorting', 'previous_feature', 'count_pass', 'count_fail'
	
	def __init__ (self, source, default_data = None, filter = (lambda x: True), verify_sorting = True):
		# look at imagesequence for an analogy
		self.source = source
		self.filter = filter
		self.default_data = default_data
		self.verify_sorting = verify_sorting
		self.previous_feature = None
		self.count_pass = 0
		self.count_fail = 0
	
	def _get_feature (self):
		'''
		this exists mainly to be replaced in subclasses
		should raise StopIteration when it reaches the end
		'''
		return next(self.source)
	
	def _try_until_pass (self):
		feature = self._get_feature()
		while not self.filter(feature):
			self.count_fail += 1
		self.count_pass += 1
		if self.default_data is None:
			return feature
		else:
			return GenomeFeature.from_genomefeature(feature, self.default_data)
	
	def __next__ (self):
		new_feature = self._try_until_pass()
		if self.verify_sorting and not (self.previous_feature is None or new_feature >= self.previous_feature): raise RuntimeError('input is not properly sorted in item %i' % self.count)
		self.previous_feature = new_feature
		return new_feature
	
	def __iter__ (self):
		return self


class SamStream (FeatureStream):
	'''
	given an iterable of pysam.AlignedSegment instances (e.g. a pysam.Samfile), yield GenomeFeatures
	self.data contains the AlignedSegment
	'''
	
	__slots__ = 'unaligned'
	
	def __init__ (self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.unaligned = 0
	
	def _get_feature (self):
		alignment = next(self.source)
		while alignment.is_unmapped:
			self.unaligned += 1
			alignment = next(self.source)
		return GenomeFeature.from_alignedsegment(alignment)


class BedStream (FeatureStream):
	'''
	given an iterable of BED-format lines (e.g. an opened BED file), yield GenomeFeatures
	you can provide a list of reference names, in order, or trust the BED data and learn automatically
	warning: if you don't provide reference names, they'll be indexed in the order they appear, so if there are any references that have no entries in the BED data the indexes won't match other data sources
	this use a lookup dictionary of reference names instead of using list.index so in theory it will perform better than GenomeFeature.from_bed
	'''
	
	__slots__ = 'fixed_references', 'reference_lookup'
	
	def __init__ (self, *args, references = None, **kwargs): # how do I pass args and kwargs while putting the new parameter last?
		super().__init__(*args, **kwargs)
		if references is None:
			self.fixed_references = False
			self.reference_lookup = collections.OrderedDict()
		else:
			self.fixed_references = True
			self.reference_lookup = collections.OrderedDict(zip(references, range(len(references))))
	
	def _get_feature (self):
		fields = next(self.source).rstrip().split()
		while len(fields) < 3: # skip bad lines
			fields = next(self.source).rstrip().split()
		reference_name, left_pos, right_pos = fields[0], int(fields[1]) + 1, int(fields[2])
		try:
			reference_id = self.reference_lookup[reference_name]
		except KeyError:
			if self.fixed_references:
				raise KeyError('unknown reference name %s in item %i' % (reference_name, self.line))
			else:
				reference_id = len(self.reference_lookup)
				self.reference_lookup[reference_name] = reference_id
		return GenomeFeature(
			reference_id =  reference_id,
			left_pos =      left_pos,
			right_pos =     right_pos,
			is_reverse =    len(fields) >= 6 and fields[5] == '-',
			data =          fields
		)
			
	@property
	def references (self):
		'''
		in case you need to check the automatically generated list
		'''
		return list(self.reference_lookup.keys())


class FastaStream (FeatureStream):
	'''
	given an iterable of FASTA-format lines (e.g. an opened FASTA file), yield GenomeFeatures of the sequences
	reference names are indexed in the order they appear, which means they're not all accessible until they're all read
	you can specify the span (number of bases per chunk) if you don't just want one at a time
	if span > 1, return_partial determines whether to return potentially shorter subsequences at the ends of input sequences
	handles buffering given the fact that FASTA is usually split across arbtirary line length
	next idea: allow overlaps!
	'''
	
	__slots__ = 'span', 'return_partial', 'references', '_sequence_buffer', 'left_pos', '_feature_generator'
	
	def __init__ (self, *args, span = 1, return_partial = True, **kwargs):
		super().__init__(*args, verify_sorting = False, **kwargs) # nonsensical to verify sorting because we're defining the coordinates
		self.span = span
		self.return_partial = return_partial
		self.references = []
		self._sequence_buffer = collections.deque()
		self.left_pos = 1 # position of the leftmost base in the sequence buffer
		self._feature_generator = self._yield_features()
	
	def _create_feature (self, length):
		seq = ''.join(self._sequence_buffer.popleft() for i in range(length))
		left_pos = self.left_pos
		self.left_pos += length
		return GenomeFeature (
			reference_id =  len(self.references) - 1, # assume we're on the most recent reference
			left_pos =      left_pos,
			right_pos =     left_pos + length - 1,
			data =          seq
		)
	
	def _purge_buffer (self):
		while len(self._sequence_buffer) >= self.span:
			yield self._create_feature(self.span)
		if len(self._sequence_buffer) > 0:
			if self.return_partial:
				yield self._create_feature(len(self._sequence_buffer))
			else:
				self._sequence_buffer.clear()			
		self.left_pos = 1
	
	def _yield_features (self):
		while True: # this seems inelegant
			while len(self._sequence_buffer) < self.span:
				try:
					new_line = next(self.source)
				except StopIteration as exception: # end of the input
					for feature in self._purge_buffer(): yield feature
					raise exception
				
				if new_line.startswith('>'): # found a reference header
					for feature in self._purge_buffer(): yield feature
					self.references += [new_line[1:].rstrip().strip()]
				else: # more sequence
					self._sequence_buffer.extend(new_line.rstrip())
			
			while len(self._sequence_buffer) >= self.span:
				yield self._create_feature(self.span)
	
	def _get_feature (self):
		return next(self._feature_generator)


class Wiggler (object):
	'''
	export GenomeFeatures to UCSC wiggle format (variableStep; https://genome.ucsc.edu/goldenPath/help/wiggle.html)
	keeps track of reference names as required
	you can add optional extra stuff for the definition line (https://genome.ucsc.edu/goldenPath/help/customTrack.html#TRACK)
	each GenomeFeature should have a numerical score as its self.data but this isn't checked (!)
	correct sorting order of the features is not checked
	make this and BED and BedGraph one class that simply sets the export format
	verify input is sorted *and* non-overlapping
	'''
	
	__slots__ = 'reference_names', 'span', 'definition', 'reference_id'
	
	def __init__ (self, reference_names, span = 1, definition = None):
		self.reference_names = reference_names
		self.span = span
		self.definition = definition
		self.reference_id = None
	
	def wiggle (self, feature):
		'''
		generates one or more lines of wiggle output from the given GenomeFeature
		if it's the first feature with a new reference_id, a new variableStep header is included
		if it's the first feature of the whole track, a track header is also included
		'''
		
		assert len(feature) == self.span
		
		result = ''
		if self.reference_id is None:
			result += 'track type=wiggle_0'
			if self.definition is not None: result += ' ' + self.definition
			result += '\n'
		if feature.reference_id != self.reference_id:
			result += 'variableStep chrom=%s' % self.reference_names[feature.reference_id]
			if self.span != 1: result += ' span=%i' % self.span
			result += '\n'
			self.reference_id = feature.reference_id
		result += '%i %s' % (feature.left_pos, feature.data)
		
		return result
		
		
class GenomePositionDeque (GenomeFeature):
	'''
	deque wrapper in which each element corresponds to a genome position
	if you add items to the ends or remove them from the ends, the coordinates change appropriately
	this can mean extending past the ends of the reference sequence!
	inserting, removing, etc. is not allowed because it's ambiguous which way to move the coordinates
	idea: it should be possible to extract a slice using only coordinates, not a GenomeFeature (with self.intersection), even though deque doesn't normally do this
	idea: don't call this a deque at all, just a sliding window or whatever
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
		or is this redundant with self.data.__getitem__ and I should really just make this the default self.__getitem__ instead?
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
	
	def count (self, value): # consider count = self.data.count but it might not work
		return self.data.count(value)
	
	def reverse (self):
		return self.data.reverse()
	
	def rotate (self, steps = 1):
		# maybe this should do popping and stuff
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
	
	def intersects(self, other):
		# this is in the wrong class
		return (
			other.reference_id == self.reference_id and
			not (other.right_pos < self.left_pos or self.right_pos < other.left_pos)
		)
	
	def intersection(self, other):
		'''
		return an instance of the overlap between this deque and another genome region, including only the data in the intersection, or None if no overlap
		use this to slice part of the original region too
		check this for bugs!
		'''
		if not self.intersects(other): return None
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
	consider making two versions: a "pull buffer" (takes as many features from the starting iterable as needed to yield one new feature at a time) and a "push buffer" (push one feature from the starting iterable at a time, yielding 0, 1, or multiple popped-off features)
	# or maybe only the "push buffer" is useful?
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

