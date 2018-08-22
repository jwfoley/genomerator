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
	def from_bed (cls, line, references):
		'''
		make a GenomeFeature from a line of a BED file
		you must give a list of reference names, in order, so it can find the index
		returns the split but unparsed fields in self.data
		consider making data an ordered dictionary for easy un-parsing
		but then also consider allowing user to specify the names of extra nonstandard fields (this will imply the number of extra fields, and then when parsing the file we will know the number of standard fields after splitting and subtracting that)
		''' 
		fields = line.rstrip().split()
		return cls(
			reference_id =  references.index(fields[0]),
			left_pos =      int(fields[1]) + 1,
			right_pos =     int(fields[2]),
			is_reverse =    len(fields) >= 6 and fields[5] == '-',
			data =          fields
		)
	
	@classmethod
	def from_gff (cls, line, references):
		'''
		make a GenomeFeature from a line of a GFF file
		you must give a list of reference names, in order, so it can find the index
		returns the split but unparsed fields in self.data
		consider making data an ordered dictionary for easy un-parsing
		'''
		fields = line.rstrip().split('\t')
		return cls(
			reference_id =  references.index(fields[0]),
			left_pos =      int(fields[3]),
			right_pos =     int(fields[4]),
			is_reverse =    fields[6] == '-',
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
	# should these be in a separate library for file formats?
	
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
		return ('%s(reference_id = %i, left_pos = %i, right_pos = %i, is_reverse = %s, data = %s)' % (self.__class__.__name__, self.reference_id, self.left_pos, self.right_pos, self.is_reverse, self.data))
