import collections
from .GenomeFeature import GenomeFeature

def read_references (reference_file):
	'''
	given a file containing reference names, in order, one per line, parse an ordered list of reference names
	if the file contains additional fields delimited by tabs, they are ignored
	'''
	references = []
	for line in reference_file:
		fields = line.rstrip().split('\t')
		if len(fields) >= 1: references.append(fields[0])
	return references

def read_reference_lengths (reference_file):
	'''
	given a file containing reference names and lengths, tab-delimited, one per line, parse an ordered dictionary of reference lengths
	if the file contains additional fields also delimited by tabs, they are ignored
	'''
	reference_lengths = collections.OrderedDict()
	for line in reference_file:
		fields = line.rstrip().split('\t')
		if len(fields) >= 2: reference_lengths[fields[0]] = int(fields[1])
	return reference_lengths

class FeatureStream (object):
	'''
	given an iterable of GenomeFeature instances, yield them back
	but also verify correct sort order, if desired, and keep a count
	you can provide a list of reference names, in order, or trust the data and learn automatically
	warning: if you don't provide reference names, they'll be indexed in the order they appear, so if there are any references that have no entries in the data the indexes won't match other data sources
	this uses a lookup dictionary of reference names instead of using list.index so in theory it will perform better than GenomeFeature's alternative constructors
	optionally replace the data with something user-specified instead
	'''
	
	__slots__ = 'source', 'default_data', 'verify_sorted', 'fixed_references', '_reference_lookup', '_feature_generator', '_previous_feature', 'count'
	
	def __init__ (self,
		source,
		references =     None,
		default_data =   None,
		verify_sorted =  False
	):
		assert isinstance(source, collections.Iterable)
		self.source = source
		self.default_data = default_data
		self.verify_sorted = verify_sorted
		self._previous_feature = None
		self.count = 0
		self._feature_generator = self._yield_features()
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
		
	def _yield_features (self):
		'''
		this exists mainly to be replaced in subclasses
		'''
		for feature in self.source: yield feature
	
	def __next__ (self):
		new_feature = next(self._feature_generator)
		if self.verify_sorted and not (self._previous_feature is None or new_feature >= self._previous_feature): raise RuntimeError('input is not properly sorted in item %i' % self.count)
		self._previous_feature = new_feature
		self.count += 1
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
	
	def __init__ (self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.unaligned = 0
	
	def _yield_features (self):
		for alignment in self.source:
			if alignment.is_unmapped:
				self.unaligned += 1
			else:
				if self.fixed_references and not alignment.reference_id == self._reference_lookup[alignment.reference_name]:
					raise RuntimeError('wrong reference ID for alignment %s' % alignment.query_name)
				else:
					yield GenomeFeature.from_alignedsegment(alignment)


class VariantStream (FeatureStream):
	'''
	given an iterable of pysam.VariantRecord instances (e.g. a pysam.VariantFile), yield GenomeFeatures
	'''
	
	def _yield_features (self):
		for variant_record in self.source:
			if self.fixed_references and not variant_record.rid == self._reference_lookup[variant_record.chrom]:
				raise RuntimeError('wrong reference ID for variant record %s' % variant_record.id)
			else:
				yield GenomeFeature.from_variantrecord(variant_record)


class BedStream (FeatureStream):
	'''
	given an iterable of BED-format lines (e.g. an opened BED file), yield GenomeFeatures
	'''
	__slots__ = 'parse'
	
	def __init__ (self, *args, parse = True, **kwargs):
		super().__init__(*args, **kwargs)
		self.parse = parse
	
	def _yield_features (self):
		for line in self.source:
			if line.startswith('#') or line.startswith('track') or line.startswith('browser'): continue
			line = line.rstrip()
			if len(line) == 0: continue
			fields = line.split()
			if len(fields) < 3: raise RuntimeError('bad format:\n%s' % line)
			yield GenomeFeature.from_bed(line = line, reference_id = self._get_reference_id(fields[0]), parse = self.parse)


class BedgraphStream (FeatureStream):
	'''
	given an iterable of BED-format lines (e.g. an opened bedGraph file), yield GenomeFeatures
	yields each genome position separately, even if file represents continguous runs, unless you specify split_spans = False
	'''
	
	__slots__ = 'split_spans'
	
	def __init__ (self, *args, split_spans = True, **kwargs):
		super().__init__(*args, **kwargs)
		self.split_spans = split_spans
	
	def _yield_features (self):
		for line in self.source:
			if line.startswith('#') or line.startswith('track') or line.startswith('browser'): continue
			line = line.rstrip()
			if len(line) == 0: continue
			fields = line.split()
			if len(fields) != 4: raise RuntimeError('bad format:\n%s' % line)
			full_feature = GenomeFeature.from_bedgraph(line = line, reference_id = self._get_reference_id(fields[0]))
			if self.split_spans:
				for position_feature in full_feature: yield position_feature
			else:
				yield full_feature
	

class WiggleStream (FeatureStream):
	'''
	given an iterable of wiggle-format lines (e.g. an opened wiggle file), yield GenomeFeatures
	yields each genome position separately, even if file represents continguous runs, unless you specify split_spans = False
	'''
	
	__slots__ = 'split_spans', '_format', '_step', '_span', '_start', '_reference_id', '_count_since_header'
	
	def __init__ (self, *args, split_spans = True, **kwargs):
		super().__init__(*args, **kwargs)
		self.split_spans = split_spans
	
	def _yield_features (self):
		for line in self.source:
			if line.startswith('#') or line.startswith('track') or line.startswith('browser'): continue
			line = line.rstrip()
			if len(line) == 0: continue
			fields = line.split()
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
			
			else: # parse a line of data
				if self._format == 'variableStep':
					if len(fields) != 2: raise RuntimeError('bad format:\n%s' % line)
					left_pos = int(fields[0])
					right_pos = left_pos + self._span - 1
					value = float(fields[1])
				if self._format == 'fixedStep':
					if len(fields) != 1: raise RuntimeError('bad format:\n%s' % line)
					left_pos = self._start + self._count_since_header * self._step
					right_pos = left_pos + self._span - 1
					value = float(fields[0])
				else:
					raise NotImplementedError('this should never happen')
				self._count_since_header += 1
				
				# split the feature into single bases if necessary
				full_feature = GenomeFeature(
					reference_id =  self._reference_id,
					left_pos =      left_pos,
					right_pos =     right_pos,
					data =          value
				)
				if self.split_spans:
					for position_feature in full_feature: yield position_feature
				else:
					yield full_feature


class GffStream (FeatureStream):
	'''
	given an iterable of GFF-format lines (e.g. an opened GFF file), yield GenomeFeatures
	optionally specify types of allowed items (e.g. gene, exon) and ignore others
	'''
	__slots__ = 'parse', 'version', 'types'
	
	def __init__ (self, *args, parse = True, version = 3, types = None, **kwargs):
		super().__init__(*args, **kwargs)
		self.parse, self.version, self.types = parse, version, types
	
	def _yield_features (self):
		for line in self.source:
			if line.startswith('#'): continue
			line = line.rstrip()
			if len(line) == 0: continue
			fields = line.split('\t')
			if len(fields) < 8: raise RuntimeError('bad format:\n%s' % line)
			result = GenomeFeature.from_gff(line = line, reference_id = self._get_reference_id(fields[0]), parse = self.parse, version = self.version)
			if self.types is None or result.data['type'] in self.types: yield result


class FastaStream (FeatureStream):
	'''
	given an iterable of FASTA-format lines (e.g. an opened FASTA file), yield GenomeFeatures of the sequences
	reference names are indexed in the order they appear, which means they're not all accessible until they're all read
	you can specify the span (number of bases per chunk) if you don't just want one at a time
	if span > 1, include_partial determines whether to return potentially shorter subsequences at the ends of input sequences
	if overlap = True, then returns overlapping sequences starting at each base position
	handles buffering given the fact that FASTA is usually split across arbtirary line length
	'''
	
	__slots__ = 'span', 'overlap', 'include_partial', '_sequence_buffer', '_left_pos', '_reference_id'
	
	def __init__ (self, *args, span = 1, overlap = False, include_partial = True, **kwargs):
		super().__init__(*args, verify_sorted = False, **kwargs) # nonsensical to verify sorting because we're defining the coordinates
		assert isinstance(span, int) and span > 0
		self.span = span
		self.overlap = overlap
		self.include_partial = include_partial
		self._sequence_buffer = collections.deque()
		self._left_pos = 1 # position of the leftmost base in the sequence buffer
		self._feature_generator = self._yield_features()
	
	def _create_feature (self, length):
		left_pos = self._left_pos
		if self.overlap:
			seq = ''.join([self._sequence_buffer.popleft()] + [self._sequence_buffer[i] for i in range(length - 1)])
			self._left_pos += 1
		else:
			seq = ''.join(self._sequence_buffer.popleft() for i in range(length))
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
			if self.include_partial:
				yield self._create_feature(len(self._sequence_buffer))
			else:
				self._sequence_buffer.clear()			
		self._left_pos = 1
	
	def _yield_features (self):
		for line in self.source:
			if line.startswith('>'): # found a reference header			
				for feature in self._purge_buffer(): yield feature
				self._reference_id = self._get_reference_id(line[1:].rstrip().strip())			
			else: # more sequence
				self._sequence_buffer.extend(line.rstrip())
				while len(self._sequence_buffer) >= self.span:
					yield self._create_feature(self.span)
		
		for feature in self._purge_buffer(): yield feature		
