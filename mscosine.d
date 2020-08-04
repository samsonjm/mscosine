/* Calculates cosine scores between MS scans.
 *
 * Author: Jonathan Samson
 * Date: 4 August 2020
 */
import std.algorithm;
import cosine;

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
	return peak_list;
}
unittest
{
	real[] mz1 = [1, 5, 3];
	real[] mz2 = [4, 2, -6];
	assert(generate_peak_list(mz1, mz2) == [-6, 1, 2, 3, 4, 5]);
	assert(generate_peak_list(mz1, mz2) == generate_peak_list(mz2, mz1));
}

int[][] create_vectors(
		int[] first_scan, 
		int[] second_scan,
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
	real[] peak_list = combine_peak_lists([1, 2, 3], [4, 5, 6]);
	//auto[real[real]] first_vector;
	//auto[real[real]] second_vector;

	//vectors = [first_vector, second_vector];
	int[][] vectors;
	//return vectors;
	return vectors;
}
unittest
{

}
