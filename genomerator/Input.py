import collections
from .GenomeFeature import GenomeFeature

class FeatureStream (object):
	'''
	given an iterable of GenomeFeature instances, yield them back
	but also verify correct sort order, if desired, and keep a count
	you can provide a list of reference names, in order, or trust the data and learn automatically
	warning: if you don't provide reference names, they'll be indexed in the order they appear, so if there are any references that have no entries in the data the indexes won't match other data sources
	this uses a lookup dictionary of reference names instead of using list.index so in theory it will perform better than GenomeFeature's alternative constructors
	optionally replace the data with something user-specified instead
	optionally filter the data to only return certain features
	'''
	
	__slots__ = 'source', 'default_data', 'filter', 'verify_order', 'fixed_references', '_reference_lookup', '_previous_feature', 'count_pass', 'count_fail'
	
	def __init__ (self,
		source,
		references =    None,
		default_data =  None,
		filter =        (lambda x: True),
		verify_order =  False
	):
		assert isinstance(source, collections.Iterable)
		self.source = source
		self.filter = filter
		self.default_data = default_data
		self.verify_order = verify_order
		self._previous_feature = None
		self.count_pass = 0
		self.count_fail = 0
		if references is None:
			self.fixed_references = False
			self._reference_lookup = collections.OrderedDict()
		else:
			self.fixed_references = True
			self._reference_lookup = collections.OrderedDict(zip(references, range(len(references))))
	
	def _get_reference_id (self, reference_name):
		'''
		look up the ID of a reference name
		if not using a predefined list, update the observed list as needed
		'''
		try:
			return self._reference_lookup[reference_name]
		except KeyError:
			if self.fixed_references:
				raise KeyError('unknown reference name %s' % (reference_name))
			else:
				reference_id = len(self._reference_lookup)
				self._reference_lookup[reference_name] = reference_id
				return reference_id
		
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
		if self.verify_order and not (self._previous_feature is None or new_feature >= self._previous_feature): raise RuntimeError('input is not properly sorted in item %i' % self.count)
		self._previous_feature = new_feature
		return new_feature
	
	def __iter__ (self):
		return self
	
	@property
	def references (self):
		'''
		in case you need to check the automatically generated list
		'''
		return list(self._reference_lookup.keys())


class SamStream (FeatureStream):
	'''
	given an iterable of pysam.AlignedSegment instances (e.g. a pysam.Samfile), yield GenomeFeatures
	'''
	
	__slots__ = 'unaligned'
	
	def _get_feature (self):
		alignment = next(self.source)
		while alignment.is_unmapped:
			self.unaligned += 1
			alignment = next(self.source)
		if self.fixed_references and not alignment.reference_id == self._reference_lookup[alignment.reference_name]:
			raise RuntimeError('wrong reference ID for alignment %s' % alignment.query_name)
		return GenomeFeature.from_alignedsegment(alignment)


class VariantStream (FeatureStream):
	'''
	given an iterable of pysam.VariantRecord instances (e.g. a pysam.VariantFile), yield GenomeFeatures
	'''
	
	def _get_feature (self):
		variant_record = next(self.source)
		if self.fixed_references and not variant_record.rid == self._reference_lookup[variant_record.chrom]:
			raise RuntimeError('wrong reference ID for variant record %s' % variant_record.id)
		return GenomeFeature.from_variantrecord(variant_record)


class BedStream (FeatureStream):
	'''
	given an iterable of BED-format lines (e.g. an opened BED file), yield GenomeFeatures
	'''
	def _get_feature (self):
		fields = next(self.source).rstrip().split()
		while fields[0].startswith('#') or len(fields) < 3: # skip bad lines
			fields = next(self.source).rstrip().split()
		return GenomeFeature(
			reference_id =  self._get_reference_id(fields[0]),
			left_pos =      int(fields[1]) + 1,
			right_pos =     int(fields[2]),
			is_reverse =    len(fields) >= 6 and fields[5] == '-',
			data =          fields
		)


class WiggleStream (FeatureStream):
	'''
	given an iterable of wiggle-format lines (e.g. an opened wiggle file), yield GenomeFeatures
	yields each genome position separately, even if file represents continguous runs, unless you specify split_spans = False
	'''
	
	__slots__ = 'split_spans', '_format', '_step', '_span', '_start', '_reference_id', '_count_since_header'
	
	def __init__ (self, *args, split_spans = True, **kwargs):
		super().__init__(*args, **kwargs)
		self.split_spans = split_spans
	
	def _get_feature (self):
		fields = next(self.source).rstrip().split()
		while len(fields) == 0 or fields[0] == 'track': fields = next(self.source).rstrip().split() # skip empty lines and track definition
		if fields[0] in ('variableStep', 'fixedStep'): # beginning of a new chromosome
			self._format = fields[0]
			assert fields[1].startswith('chrom=')
			self._reference_id = self._get_reference_id(fields[1][6:])
			
			# parse format-specific settings
			if self._format == 'variableStep':
				if len(fields) > 2:
					assert len(fields) == 3
					assert fields[2].startswith('span=')
					self._span = int(fields[2][5:])
				else:
					self._span = 1
			elif self._format == 'fixedStep':
				assert fields[2].startswith('start=')
				self._start = int(fields[2][6:])
				assert fields[3].startswith('step=')
				self._step = int(fields[3][5:])
				if len(fields) > 4:
					assert len(fields) == 5
					assert fields[4].startswith('span=')
					self._span = int(fields[4][5:])
				else:
					self._span = 1				
			else:
				raise NotImplementedError('this should never happen')
			
			self._count_since_header = 0
			fields = next(self.source).rstrip().split() # read another line, which we assume (!) is now data
		
		# parse a line of data
		if self._format == 'variableStep':
			assert len(fields) == 2
			left_pos = int(fields[0])
			right_pos = left_pos + self._span - 1
			value = float(fields[1])
		elif self._format == 'fixedStep':
			assert len(fields) == 1
			left_pos = self._start + self._count_since_header * self._step
			right_pos = left_pos + self._span - 1
			value = float(fields[0])
		else:
			raise NotImplementedError('this should never happen')
		
		self._count_since_header += 1
		return GenomeFeature(
			reference_id =  self._reference_id,
			left_pos =      left_pos,
			right_pos =     right_pos,
			data =          value
		)



class GffStream (FeatureStream):
	'''
	given an iterable of GFF-format lines (e.g. an opened GFF file), yield GenomeFeatures
	'''
	
	def _get_feature (self):
		fields = next(self.source).rstrip().split('\t')
		while fields[0].startswith('#') or len(fields) < 8: # skip bad lines
			fields = next(self.source).rstrip().split('\t')
		return GenomeFeature(
			reference_id =  self._get_reference_id(fields[0]),
			left_pos =      int(fields[3]),
			right_pos =     int(fields[4]),
			is_reverse =    fields[6] == '-',
			data =          fields
		)


class FastaStream (FeatureStream):
	'''
	given an iterable of FASTA-format lines (e.g. an opened FASTA file), yield GenomeFeatures of the sequences
	reference names are indexed in the order they appear, which means they're not all accessible until they're all read
	you can specify the span (number of bases per chunk) if you don't just want one at a time
	if span > 1, return_partial determines whether to return potentially shorter subsequences at the ends of input sequences
	handles buffering given the fact that FASTA is usually split across arbtirary line length
	'''
	
	__slots__ = 'span', 'return_partial', '_sequence_buffer', '_left_pos', '_reference_id', '_feature_generator'
	
	def __init__ (self, *args, span = 1, return_partial = True, **kwargs):
		super().__init__(*args, verify_order = False, **kwargs) # nonsensical to verify sorting because we're defining the coordinates
		self.span = span
		self.return_partial = return_partial
		self._sequence_buffer = collections.deque()
		self.left_pos = 1 # position of the leftmost base in the sequence buffer
		self._feature_generator = self._yield_features()
	
	def _create_feature (self, length):
		seq = ''.join(self._sequence_buffer.popleft() for i in range(length))
		left_pos = self._left_pos
		self._left_pos += length
		return GenomeFeature (
			reference_id =  self._reference_id,
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
					line = next(self.source)
				except StopIteration as exception: # end of the input
					for feature in self._purge_buffer(): yield feature
					raise exception
				
				if line.startswith('>'): # found a reference header
					for feature in self._purge_buffer(): yield feature
					self._reference_id = self._get_reference_id(line[1:].rstrip().strip())
				else: # more sequence
					self._sequence_buffer.extend(line.rstrip())
			
			while len(self._sequence_buffer) >= self.span:
				yield self._create_feature(self.span)
	
	def _get_feature (self):
		return next(self._feature_generator)
