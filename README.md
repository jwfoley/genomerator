# genomerator

Python generators for iterating over genome data

## Glossary

* **Reference**: A single physical reference sequence for a genome, transcriptome, gene, etc. with locations described by integer coordinates (e.g. base positions). Each chromosome or contig is one reference. A genome may contain many references. For ease of data management, references are tracked by **reference ID (index)**, a non-negative integer expressing each reference's rank in your ordering system. The first reference's index is 0, because these indexes are used to access the list of reference names and Python starts list numbering at 0. But the order is arbitrary; you can put `chr10` before `chr2` for alphabetical sorting or after it for numerical consistency, as long as you always keep the references in the same order.
* **Feature**: Any sort of annotation with coordinates on a reference. Genes, exons, SNPs, enriched regions, etc. are all kinds of features. In this module, a feature can only be paired with one reference.
* **Left and right positions**: The "left" position of a feature is the base position within it that has the lowest coordinate (the farthest left as you look at it on a genome browser), expressed as a positive integer, and "right" is the highest coordinate. Other software often refers to these as "start" and "stop"/"end" positions, but that is ambiguous depending on the orientation of the feature. Note: in this module, unlike some other software, both left and right positions are always 1-based, i.e. the first base of the reference is base #1. Thus if the feature is only one base long, like a SNP, its left and right positions are the same number. This means genome coordinates inside the module will match coordinates on the genome browser.
* **Strand orientation**: A feature may be oriented in the same direction as the reference (forward strand, also called plus strand or Watson strand) or the opposite direction (reverse strand, minus strand, Crick strand). E.g. if a gene is on the forward strand, its transcription initiation site has a lower-numbered position than its transcription termination site. Some kinds of features do not have any orientation; for simplicity these are treated as forward-oriented.
* **Start and end positions**: In this module, the "start" and "end" of a feature are designated relative to the orientation of the feature, not the orientation of the reference. So if the feature is forward-oriented, the start position is the left position and the end position is the right position. But if it is reverse-oriented, the start position is the right position and the end position is the left position. It is important to keep these terms distinct as each context is useful in different situations.



## `GenomeFeature`

This is a general-purpose container for any kind of data with reference coordinates. It always contains several important attributes:

* `reference_id`: The index number of the reference this feature is on (first reference is number 0).
* `left_pos`, `right_pos`: The leftmost and rightmost position numbers of the feature on that reference (first base of the reference is always position 1 regardless of whether it's the left or right end of a feature).
* `is_reverse`: A Boolean value describing the orientation of the feature relative to the reference. If `True`, it is reverse-oriented. If `False`, it is forward-oriented, or has no orientation.
* `data`: This attribute can contain any kind of information regarding the feature: a count, a score, a name, even a large data structure containing an unlimited amount of information.

### Constructor: `GenomeFeature(reference_id, left_pos, right_pos = None, is_reverse = False, data = None`

The `reference_id` and `left_pos` must always be specified. If no `right_pos` is specified, it is assumed to be the same as the `left_pos` for a single-base feature, e.g. a SNP.

Several alternative constructors exist for creating a `GenomeFeature` from certain kinds of data:

* `GenomeFeature.from_genomefeature(other, data = None)`: Create a new instance at the same location as the old one, but you can specify new data for it to contain (if you don't specify, it has the same data as the old instance).
* `GenomeFeature.from_alignedsegment(aligned_segment)` and `GenomeFeature.from_variantrecord(variant_record)`: These convert a `pysam.AlignedSegment` or `pysam.VariantRecord`, respectively, into a `GenomeFeature`. The original object is stored in the `data` attribute. E.g. you can still access the data of a converted `AlignedSegment` as `my_feature.data.sequence`.
* `GenomeFeature.from_bed(line, reference_id)` and the `.from_bedgraph`, `.from_gff` methods with the same syntax: These create a new `GenomeFeature` from a single line of a UCSC BED, UCSC bedGraph, or GFF file, respectively. `.from_bedgraph` puts the score (as a `float`) in `data`. `from_bed` and `from_gff` have another argument, `parse` (default true): if enabled, all the data fields defined in the format are parsed by name into a dictionary stored in `data`, or if disabled to save computing time, all the delimiter-split fields (as a list of strings) are stored there instead.

### Extracting positions

`GenomeFeature.start_pos` and `GenomeFeature.end_pos` each return an `int` of the feature's start or end position, respectively, according to its strand orientation.

`GenomeFeature.left()`, `.right()`, `.start()`, and `.end()` each return a new, complete `GenomeFeature` instance consisting only of the leftmost/rightmost/start/end position of the original instance, but containing the same data.

`my_feature[index]` returns a new instance with the same data at the specified position, where `index` is relative to the left position of the original instance. For example, if `my_feature`'s left position is 2355 and right position is 2358, `my_feature[0]` returns a new instance with left 2355, right 2355; `my_feature[2]` is left 2357, right 2357. Reverse indexes and contiguous forward-oriented slices are also allowed: `my_feature[-1]` is left 2358, right 2358, and `my_feature[2:]` is left 2357, right 2358. `GenomeFeature.get_pos(position)` does the same thing but with a genome position rather than a relative offset index, e.g. `my_feature[2]` is the same as `my_feature.get_pos(2357)`. But note that, consistent with Python syntax, the stop position of a slice is still one past the end of the desired section; so `my_feature.get_pos(2355:2359)` will return the entire feature, length 4.

`len(my_feature)` returns the number of positions covered by a feature, e.g. `my_feature` has a length of 4. A single-position feature, e.g. left 2355 and right 2355, has a length of 1.

Iterating over a `GenomeFeature` yields single-base features from left to right, e.g. `for i in my_feature:` gives you position 2355, then 2356, then 2357, then 2358.

### Modifying positions

In addition to simply setting `left_pos` and `right_pos`, you can also perform operations to shift one or both of them a specific distance. `GenomeFeature.shift_left(distance)` shifts only the left position by the specified distance and `.shift_right` shifts only the right position. Note that both methods shift the feature to the right (adding `distance` to `left_pos` or `right_pos`), unless `distance` is negative, and then they shift to the left. You cannot shift them so far that they become backwards and `right_pos` is less than `left_pos`. `GenomeFeature.shift(distance)` shifts both coordinates by the same distance. No amount of shifting will change the `reference_id`; you can shift past the right end of the reference, with position numbers greater than its entire length, or past the left end, with negative positions, if you have some reason to do that.

There are also versions that shift relative to the orientation of the feature rather than the orientation of the reference. `shift_start` and `shift_end` shift the start and end positions, but note that they shift in the direction of the feature. For example, if the feature is forward-oriented, `shift_start` moves its *left* position to the *right* (like `shift_left`); but if it is reverse-oriented, `shift_start` moves its *right* position to the *left*. `shift_forward` shifts both coordinates by the same distance, still in the direction of the feature.

For convenience, the `+` and `-` operators essentially apply `shift`, but create a new instance, to keep syntax consistent. So `GenomeFeature(reference_id = 5, left_pos = 23, right_pos = 28) + 2` gives you a new feature with left position 25, right 30.

`GenomeFeature.switch_strand()` simply changes `is_reverse` to the opposite of what it was before.

### Comparing features

The `<`, `>`, and `==` operators compare `GenomeFeature` instances' positions **by sorting order**. First they compare by `reference_id`: all features on reference 2 are "less than" features on reference 3. If the features are on the same reference, they are compared **only by left position**: a feature with left position 23, right position 42 is less than a feature with left 25, right 38. Strand orientation is ignored by these operators. This is arbitrary and somewhat unintuitive, but it allows features to be put in conventional sorting order with standard Python syntax: the build-in `sorted` function will sort a list of features in the same order expected by `samtools` and virtually all other genomics software.

More intuitively, you may want to check whether one feature is entirely to the left of another feature, with no overlap. `GenomeFeature.left_of(other)` and `GenomeFeature.right_of(other)` do this. The former will return `True` if the given feature's right position is to the left of the other feature's left position (or this feature is on an earlier reference than the other).

`GenomeFeature.same_as(other)` compares whether two features have the same reference, left position, right position, and orientation. They may still contain different data.

`GenomeFeature.intersects(other)` tests whether two features overlap at any position. They must be on the same reference to overlap. If they do overlap, `GenomeFeature.intersection(other)` returns a new instance containing only the overlap region, and the same data as the feature it was called from (the data from `other` are discarded). Contrariwise, `GenomeFeature.union(other)` returns a new instances containing the longest possible region covered by either feature, as long as there is no gap between then.

### Progress calculation

When you are iterating over a whole genome, you may wish to check how far you've gotten, e.g. to estimate remaining time until your script finishes running. `GenomeFeature.progress(reference_lengths)` calculates the proportion of the whole genome that lies before the given feature, if you give it a list of the lengths of all references.

### Exporting

The `GenomeFeature` class contains a few methods to export an instance into a single line of a common genomics file format.

`GenomeFeature.bed(reference_names, name = None, score = None, thick_range = None, rgb = None, block_count = None, block_sizes = None, block_starts = None, use_strand = True)` exports into UCSC BED format. You must provide the list of reference names, in the correct order, so it can convert `reference_id` to a name. You can provide data for the optional fields if desired, otherwise they are filled in with `.` per the specification. If `use_strand` is disabled, the strand field will also be filled with `.` to indicate a nondirectional feature.

`GenomeFeature.gff(reference_names, type, source = None, score = None, phase = None, attributes = None, use_strand = True)` behaves similarly for GFF format, but one of the fields, `type`, is required by the format specification so this method requires it too. `attributes` should be a string containing the entire list of attributes for that field in the format.

`GenomeFeature.bedgraph(reference_names)` exports to the simpler UCSC bedGraph format, so it only needs the list of reference names. It assumes the feature's `data` attribute contains the numerical score for this genome region.


## Specialized subclasses

### `GenomeArray`

This subclass stores position-specific data for a genome region. That is, its `data` attribute is an iterable containing some kind of data value for every position in the range. Thus `len(my_genomearray)` is always the same as `len(my_genomearray.data)`. Like `collections.defaultdict`, `GenomeArray` is instantiated with a special argument, `default_factory`, which is a function that creates the default data value for positions in the array. The default `default_factory` is `float`, which will create a data value of 0.0 for every position. You can change it to e.g. `list` or `dict` and then each position will have its own separate container for larger-scale data.

`GenomeArray`'s methods for accessing individual positions within the region are similar to base `GenomeFeature`'s, except they return instances containing position-specific data. For example, if `my_genomearray` spans from left position 13 to right position 20, `my_genomearray.data` has length 8, and `my_genomearray[:3]` returns a new instance containing only positions 13 through 15 and the data values for only those first 3 positions. 

Likewise the methods for changing coordinates also change the data. `shift_left`, `shift_right`, `shift_start`, and `shift_end` either delete data from the corresponding end of the array (if the shift shortens it) or add new values from the default factory (if lengthening it). `shift` and `shift_forward` remove from one end and add to the other, so the same data values will still correspond to the same positions. If you don't want to change the data while moving the array, `move(distance)` is like `shift(distance)` but keeps the data intact. 

If you want to add specific data while changing the coordinates, the methods are analogous to `collections.deque`: `append(value)`, `extend(values)`, and `pop()` operate on the right end and `appendleft(value)`, `extendleft(values)`, and `popleft(values)` operate on the left end. There is even `rotate(steps = 1)`, which moves the array while simply rotating the contents of the data, so the values from the end that's being shortened are added to the end being lengthened. `reverse()` reverses the order of the data. `intersection(other)` returns only the data values from the intersection region.

Because the data value is an iterable, you can use `count(value)` directly on the array to count the number of positions whose data is `value`.

