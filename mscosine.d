/* Calculates cosine scores between MS scans.
 *
 * Author: Jonathan Samson
 * Date: 4 August 2020
 */
import std.algorithm;
import cosine;
import scans;

real[] combine_peak_lists(real[] mz1, real[] mz2)
/* Creates a combined list of peaks from the separate peak lists.
 *
 * Arguments:
 *	mz1 - The list of mass/charge ratios from one scan.
 *	mz2 - The list of mass/charge ratios from another scan.
 * Returns:
 *	peak_list - The sorted combined list of mass/charge ratios.
 */
{
	real[] peak_list = mz1;
	foreach(peak; mz2)
		peak_list ~= peak;
	peak_list.sort();
	peak_list.length -= peak_list.uniq().copy(peak_list).length;
	return peak_list;
}
unittest
{
	real[] mz1 = [1, 5, 3];
	real[] mz2 = [3, 2, -6];
	assert(combine_peak_lists(mz1, mz2) == [-6, 1, 2, 3, 5]);
	assert(combine_peak_lists(mz1, mz2) == combine_peak_lists(mz2, mz1));
}

real[real][2] create_vectors(
		real[real] first_scan, 
		real[real] second_scan,
		real threshold = 0.01)
/* Creates vectors of equal length for each scan.
 *
 * Arguments:
 *	first_scan - scans with mz:intensity.
 *	second_scan - scans with mz:intensity.
 *	threshold - mz tolerance to be considered the same peak.
 * Returns:
 *	vectors - first scan then second scan vectors.
 *
 * If the scan has a peak, a tuple will be present with mz:intensity.
 * If more than 1 peak are within the threshold within the same scan, 
 * the peak intensities will be added.  If no peak has been found in 
 * the scan, the intensity will be 0.
 *
 * The algorithm works from low m/z to high in a greedy fashion; the
 * threshold comparison will be made to the smallest peak.
 */
{
	real[] first_peaks = first_scan.keys();
	real[] second_peaks = second_scan.keys();
	real[] peak_list = combine_peak_lists(first_peaks, second_peaks);

	real[real] first_vector;
	real[real] second_vector;
	sort(first_peaks);
	sort(second_peaks);
	real[] to_pop;
	while(peak_list.length > 0)
	{
		while(first_peaks.length > 0 && 
		      peak_list[0] + threshold >= first_peaks[0])
		{
			to_pop ~= first_peaks[0];
			first_vector[peak_list[0]] = first_scan[first_peaks[0]];
			first_peaks = first_peaks[1..$];
		}
		while(second_peaks.length > 0 &&
		      peak_list[0] + threshold >= second_peaks[0])
		{
			to_pop ~= second_peaks[0];
			second_vector[peak_list[0]] = second_scan[second_peaks[0]];
			second_peaks = second_peaks[1..$];
		}
		to_pop.length -= to_pop.uniq().copy(to_pop).length;
		peak_list = peak_list[(to_pop.length)..$];
		to_pop = [];
	}
	real[real][2] vectors = [first_vector, second_vector];
	return vectors;
}
unittest
{
	real[real] scan1_peaks = [
		100.1: 10000.100,
		200.1: 15000.2,
		300.1: 16000.3
	];
	real[real] scan2_peaks = [
		100.1: 10000.100,
		200.1: 15000.2,
		300.1: 16000.3
	];
	assert(create_vectors(scan1_peaks, scan2_peaks) == 
			[scan1_peaks, scan2_peaks]);
}
