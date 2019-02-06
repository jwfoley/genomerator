#! /usr/bin/env python3

import sys, click, genomerator

@click.command()
@click.argument('reference_list', type = click.File('r'))
@click.argument('in_bed', required = False, type = click.File('r'), default = click.open_file('-', 'r'))
@click.argument('out_bed', required = False, type = click.File('w'), default = click.open_file('-', 'w'))
def sort_bed (reference_list, in_bed, out_bed):
	'''
	given a list of reference names, sort a BED file into the same order
	features on unrecognized references are discarded
	'''
	
	feature_stream = genomerator.BedStream(
		in_bed,
		assert_sorted = False,
		parse = False,
		references = genomerator.read_references(reference_list)
	)
	for feature in sorted(feature_stream): out_bed.write('\t'.join(map(str, feature.data)) + '\n')
	if feature_stream.count_discarded > 0: print('warning: %i features discarded because reference wasn\'t in the list' % feature_stream.count_discarded, file = sys.stderr)	
	
if __name__ == '__main__':
	sort_bed()

