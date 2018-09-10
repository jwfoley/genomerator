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
		result += '%i %s\n' % (feature.left_pos, feature.data)
		
		return result
