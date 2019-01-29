import click

class ProgressBar (object):
	'''
	wrapper for a click.progressbar that keeps track of reference lengths and updates with GenomeFeatures
	__del__ isn't guaranteed to be called at the end of a script so it helps to explicitly .finish when done
	'''
	
	__slots__ = '_progressbar', 'reference_lengths', 'progress'
	
	def __init__ (self, reference_lengths):
		self._progressbar = click.progressbar(length = 1)
		self.reference_lengths = reference_lengths
		self.progress = 0
	
	def update (self, feature):
		new_progress = feature.progress(self.reference_lengths)
		self._progressbar.update(new_progress - self.progress)
		self.progress = new_progress
	
	def finish (self):
		self._progressbar.update(1)
		self._progressbar.render_finish()
	
	def __del__ (self):
		self._progressbar.render_finish()

