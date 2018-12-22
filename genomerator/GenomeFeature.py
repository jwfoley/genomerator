import copy

class GenomeFeature (object):
	'''
	container for a generic genome feature with reference ID number (not name), left and right positions (1-based), strand (as 'is_reverse'), and embedded data of any kind
	if no right_pos is specified, it represents a single genome position (and right_pos = left_pos)
	non-stranded features default to the forward strand (is_reverse = False)
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
		wrap a pysam.AlignedSegment into a GenomeFeature
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
	def from_variantrecord (cls, variant_record):
		'''
		wrap a pysam.VariantRecord into a GenomeFeature
		the original VariantRecord is still stored as self.data
		'''
		fields = line.rstrip().split('\t')
		return cls(
			reference_id =  variant_record.rid,
			left_pos =      variant_record.start + 1,
			right_pos =     variant_record.stop,
			data =          variant_record
		)
	
	@classmethod
	def from_bed (cls, line, reference_id, parse = True):
		'''
		make a GenomeFeature from a line of a BED file
		https://genome.ucsc.edu/FAQ/FAQformat.html#format1
		you must give a reference_id because this format only has the name
		parses the fields into self.data if parse = True
		otherwise simply puts the unparsed, split strings into self.data (faster)
		''' 
		fields = line.rstrip().split()
		if parse:
			data = {}
			try:
				data['name'] =         None if fields[3] == '.' else fields[3]
				data['score'] =        None if fields[4] == '.' else float(fields[4])
				data['thickStart'] =   None if fields[6] == '.' else int(fields[6]) + 1
				data['thickEnd'] =     None if fields[7] == '.' else int(fields[7])
				data['itemRgb'] =      None if fields[8] == '.' else list(int(x) for x in fields[8].split(','))
				data['blockCount'] =   None if fields[9] == '.' else int(fields[9])
				data['blockSizes'] =   None if fields[10] == '.' else list(int(x) for x in fields[10].split(','))
				data['blockStarts'] =  None if fields[11] == '.' else list(int(x) for x in fields[11].split(','))
			except IndexError: # keep parsing until we run out of columns, then just stop with what we have
				pass
		else:
			data = fields
		return cls(
			reference_id =  reference_id,
			left_pos =      int(fields[1]) + 1,
			right_pos =     int(fields[2]),
			is_reverse =    len(fields) >= 6 and fields[5] == '-',
			data =          data
		)
	
	@classmethod
	def from_bedgraph (cls, line, reference_id):
		'''
		make a GenomeFeature from a line of a bedGraph file
		https://genome.ucsc.edu/FAQ/FAQformat.html#format1.8
		you must give a reference_id because this format only has the name
		parses the dataValue into self.data
		'''
		fields = line.rstrip().split()
		return cls(
			reference_id =  reference_id,
			left_pos =      int(fields[1]) + 1,
			right_pos =     int(fields[2]),
			data =          float(fields[3])
		)
	
	@classmethod
	def from_gff (cls, line, reference_id, parse = True, version = 3):
		'''
		make a GenomeFeature from a line of a GFF file
		https://genome.ucsc.edu/FAQ/FAQformat.html#format3
		you must give a list of reference names, in order, so it can find the index
		parses the fields into self.data if parse = True
		otherwise simply puts the unparsed, split strings into self.data (faster)
		'''
		fields = line.rstrip().split('\t')
		if parse:
			data = {
				'source':      None if fields[1] == '.' else fields[1],
				'type':        fields[2],
				'score':       None if fields[5] == '.' else float(fields[5]),
				'phase':       None if fields[7] == '.' else int(fields[7])
			}
			for attribute in fields[8].split(';'):
				if len(attribute) == 0: continue # some people leave trailing semicolons
				if version == 2:
					attribute = attribute.strip().rstrip() # possible leading/trailing spaces if there are spaces between items
					first_space_pos = attribute.index(' ') # can't simply split because the value may also contain spaces
					tag, value = attribute[:first_space_pos], attribute[first_space_pos + 1:]
					if value.startswith('"') and value.endswith('"'): value = value[1:-1]
				elif version == 3:
					tag, value = attribute.split('=')
				else:
					raise RuntimeError('unknown GFF version %s' % str(version))
				data[tag] = value
		else:
			data = fields
		return cls(
			reference_id =  reference_id,
			left_pos =      int(fields[3]),
			right_pos =     int(fields[4]),
			is_reverse =    fields[6] == '-',
			data =          data
		)
	
	# accessing the attributes
	
	@property
	def left (self):
		'''
		return a new instance containing only the leftmost position, for quick comparisons
		'''
		return GenomeFeature(reference_id = self.reference_id, left_pos = self.left_pos, is_reverse = self.is_reverse, data = self.data)
	
	@property
	def right (self):
		'''
		return a new instance containing only the rightmost position, for quick comparisons
		'''
		return GenomeFeature(reference_id = self.reference_id, left_pos = self.right_pos, is_reverse = self.is_reverse, data = self.data)
	
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
	
	def __len__ (self):
		return self.right_pos - self.left_pos + 1
	

	# extracting positions
	
	def _compute_position (self, index):
		'''
		helper function to convert an offset index to a base position
		position can be between the left position and 1 past the right position (for slicing)
		'''
		if 0 <= index <= len(self):
			return self.left_pos + index
		elif -1 >= index >= -len(self):
			return self.right_pos + index + 1
		else:
			raise IndexError
	
	def __getitem__ (self, index):
		'''
		return a new instance containing the same data but only at the specified position (as an offset index or slice)
		'''
		if isinstance(index, int): # single index
			return GenomeFeature(reference_id = self.reference_id, left_pos = self._compute_position(index), is_reverse = self.is_reverse, data = self.data)
				
		elif isinstance(index, slice): # slice
			if index.step is not None: raise IndexError('can only slice by step = 1')
			if index.start is not None and index.stop is not None and index.stop <= index.start: raise IndexError
			
			new_left = (self.left_pos if index.start is None else self._compute_position(index.start))
			new_right = (self.right_pos if index.stop is None else self._compute_position(index.stop) - (index.stop >= 0))
			
			return GenomeFeature(reference_id = self.reference_id, left_pos = new_left, right_pos = new_right, is_reverse = self.is_reverse, data = self.data)

		else:
			raise TypeError
	
	def get_pos (self, position):
		'''
		return a new instance containing the same data but only at the specified position (as a genome position)
		'''
		if isinstance(position, int):
			if not self.left_pos <= position <= self.right_pos: raise IndexError
			return self[position - self.left_pos]
		elif isinstance(position, slice):
			start = 0 if position.start is None else position.start - self.left_pos
			if not 0 <= start < len(self): raise IndexError
			stop = len(self) if position.stop is None else position.stop - self.left_pos
			if not 0 < stop <= len(self): raise IndexError
			return self[start:stop]
		else:
			raise TypeError

	
	def __iter__ (self):
		'''
		return a new instance for each single genome position
		'''
		for index in range(len(self)): yield self[index]
		
	
	# modifying the coordinates
	
	def shift_left (self, distance):
		'''
		shift the left coordinate by the specified distance (note: positive values shift it to the right)
		'''
		if distance > len(self) - 1: raise IndexError
		self.left_pos += distance
	
	def shift_right (self, distance):
		'''
		shift the right coordinate by the specified distance
		'''
		if distance < -(len(self) - 1): raise IndexError
		self.right_pos += distance
	
	def shift (self, distance):
		'''
		shift both coordinates by the specified distance
		'''
		self.left_pos += distance
		self.right_pos += distance
	
	def shift_start (self, distance):
		'''
		shift the beginning of the feature by the specified distance, strand-specifically
		so if self.is_reverse is false, this shifts the left position to the right
		or if self.is_reverse is true, this shifts the right position to the left				
		'''
		if self.is_reverse:
			self.shift_right(-distance)
		else:
			self.shift_left(distance)
	
	def shift_end (self, distance):
		'''
		analogous to shift_start
		'''
		if self.is_reverse:
			self.shift_left(-distance)
		else:
			self.shift_right(distance)
	
	def shift_forward (self, distance):
		'''
		shift both coordinates by the specified distance, strand-specifically
		so is self.is_reverse is false, it moves right, otherwise it moves left
		unless distance is negative, in which case it's the opposite
		'''
		if self.is_reverse:
			self.shift(-distance)
		else:
			self.shift(distance)	

	def switch_strand (self):
		self.is_reverse = not self.is_reverse	

	
	# returning a copy with modified coordinates	
	
	def shifted_left (self, distance):
		new_instance = copy.copy(self)
		new_instance.shift_left(distance)
		return new_instance
	
	def shifted_right (self, distance):
		new_instance = copy.copy(self)
		new_instance.shift_right(distance)
		return new_instance
	
	def shifted (self, distance):
		new_instance = copy.copy(self)
		new_instance.shift(distance)
		return new_instance
		
	def shifted_start (self, distance):
		new_instance = copy.copy(self)
		new_instance.shift_start(distance)
		return new_instance
		
	def shifted_end (self, distance):
		new_instance = copy.copy(self)
		new_instance.shift_end(distance)
		return new_instance
	
	def shifted_forward (self, distance):
		new_instance = copy.copy(self)
		new_instance.shift_forward(distance)
		return new_instance
	
	def __add__ (self, distance):
		'''
		same as shifted
		'''
		return self.shifted(distance)

	def __sub__ (self, distance):
		'''
		opposite of __add__
		'''
		return self.__add__(-distance)
	
	def switched_strand (self):
		new_instance = copy.copy(self)
		new_instance.switch_strand()
		return new_instance
		
		
	# logical comparisons
	
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
	
	def intersects (self, other):
		'''
		test whether two features overlap
		'''
		return other.reference_id == self.reference_id and not (other.right_pos < self.left_pos or self.right_pos < other.left_pos)
	
	
	# distance comparisons
	# these return infinity if features aren't on the same reference
	# positive infinity if other feature is on a higher-number reference, negative infinity if on a lower-numbered one
	
	def left_offset (self, other):
		'''
		compute offset of other feature's left position relative to this feature's left position
		sign of result corresponds to reference orientation
		'''
		if other.reference_id > self.reference_id:
			return float('inf')
		elif other.reference_id < self.reference_id:
			return float('-inf')
		else:
			return other.left_pos - self.left_pos

	def right_offset (self, other):
		'''
		compute offset of other feature's right position relative to this feature's right position
		sign of result corresponds to reference orientation
		'''
		if other.reference_id > self.reference_id:
			return float('inf')
		elif other.reference_id < self.reference_id:
			return float('-inf')
		else:
			return other.right_pos - self.right_pos
	
	def start_offset (self, other):
		'''
		compute offset of other feature's start position relative to this feature's start position
		sign of result corresponds to this feature's orientation, not reference orientation
		'''
		if other.reference_id > self.reference_id:
			return float('inf')
		elif other.reference_id < self.reference_id:
			return float('-inf')
		elif self.is_reverse:
			return self.start_pos - other.start_pos
		else:
			return other.start_pos - self.start_pos
	
	def distance_to (self, other):
		'''
		compute offset of other feature's nearest position relative to this feature's nearest position
		if they intersect, distance is zero
		sign of result corresponds to reference orientation
		'''
		if other.reference_id > self.reference_id:
			return float('inf')
		elif other.reference_id < self.reference_id:
			return float('-inf')
		elif self.intersects(other):
			return 0
		elif other.left_of(self):
			return other.right_pos - self.left_pos
		else:
			return other.left_pos - self.right_pos
	
	def __sub__(other, self):
		'''
		compute distance from one feature to another
		'''
		return self.distance_to(other)
	
	
	# set-like operations
	
	def intersection (self, other):
		'''
		return a GenomeFeature of the overlap between two regions (preserves this one's data), or None if no overlap
		'''
		if not self.intersects(other): return None
		return GenomeFeature(
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
	
	def space_between (self, other):
		'''
		return a GenomeFeature of the region between two GenomeFeatures
		returns None if they are on different references or they have no gap between them
		'''
		if (
			other.reference_id != self.reference_id or
			self.intersects(other)
		):
			return None
		else:
			if self.left_of(other):
				if self.right_pos + 1 == other.left_pos: return None
				left_pos = self.right_pos + 1
				right_pos = other.left_pos - 1
			else:
				if other.right_pos + 1 == self.left_pos: return None
				left_pos = other.right_pos + 1
				right_pos = self.left_pos - 1
			return GenomeFeature (
				reference_id =  self.reference_id,
				left_pos =      left_pos,
				right_pos =     right_pos,
				data =          self.data
			)
	
	# progress (proportion of genome traversed)
	
	def progress (self, reference_lengths):
		total_length = sum(reference_lengths)
		assert total_length > 0
		traversed_length = sum(reference_lengths[:self.reference_id]) + self.left_pos
		return traversed_length / total_length
	
	
	# exporting
	
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
	
	def gff (self,
		reference_names,
		type,
		source =      None,
		score =       None,
		phase =       None,
		attributes =  None,
		use_strand =  True
	):
		'''
		format the feature as a GFF3 line
		requires a list of reference names since the feature only knows its index
		type is required by the GFF spec but other fields are optional
		'''
		assert type != 'CDS' or phase is not None # phase is required for all CDS features
		return '\t'.join([
			reference_names[self.reference_id],                        # seqid
			source if source is not None else '.',                     # source
			type,                                                      # type
			str(self.left_pos),                                        # start
			str(self.right_pos),                                       # end
			str(score) if score is not None else '.',                  # score
			('-' if self.is_reverse else '+') if use_strand else '.',  # strand
			str(phase) if phase is not None else '.',                  # phase
			attributes if attributes is not None else '.'              # attributes
		])
		
	def __repr__ (self):
		return ('%s(reference_id = %i, left_pos = %i, right_pos = %i, is_reverse = %s, data = %s)' % (self.__class__.__name__, self.reference_id, self.left_pos, self.right_pos, self.is_reverse, self.data.__repr__()))
