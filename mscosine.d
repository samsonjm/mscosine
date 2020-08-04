/* Calculates cosine scores between MS scans.
 *
 * Author: Jonathan Samson
 * Date: 4 August 2020
 */
module mscosine;
import std.algorithm;
import cosine;
import std.math;
import std.exception;

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
		real total_value = 0;
		while(first_peaks.length > 0 && 
		      peak_list[0] + threshold >= first_peaks[0])
		{
			to_pop ~= first_peaks[0];
			total_value += first_scan[first_peaks[0]];
			first_peaks = first_peaks[1..$];
		}
		first_vector[peak_list[0]] = total_value;
		total_value = 0;
		while(second_peaks.length > 0 &&
		      peak_list[0] + threshold >= second_peaks[0])
		{
			to_pop ~= second_peaks[0];
			total_value += second_scan[second_peaks[0]];
			second_peaks = second_peaks[1..$];
		}
		second_vector[peak_list[0]] = total_value;
		sort(to_pop);
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
		100.2: 500,
		200.1: 15000.2,
		300.1: 16000.3
	];
	real[real] scan2_peaks = [
		100.1: 10000.1,
		200.2: 15000.2,
		300.3: 16000.3
	];
	real[real] scan1_vector = [
		100.1: 10500.100,
		200.1: 15000.2,
		300.1: 16000.3,
		300.3: 0 
	];
	real[real] scan2_vector = [
		100.1: 10000.1,
		200.1: 15000.2,
		300.1: 0,
		300.3: 16000.3
	];
	assert(create_vectors(scan1_peaks, scan2_peaks, 0.1) == [
		scan1_vector, scan2_vector]);
}

bool not_all_zeroes(real[] list)
/* Checks whether all values are 0.
 * Arguments:
 *	list - A list of real values.
 * Returns:
 *	non_zeroes - true only if there is at least 1 non-zero value.
 */
{
	bool non_zeroes = false;
	foreach(real value; list)
	{
		if(value != 0)
		{
			non_zeroes = true;
			break;
		}
	}
	return non_zeroes;
}
unittest
{
	real[] test = [0, 0.0, 0.000000];
	assert(not_all_zeroes(test) == false);
	test = [0, 0.0, 0.00000001];
	assert(not_all_zeroes(test) == true);
}

real find_cosine_score(
		real[real] first_scan, 
		real[real] second_scan, 
		real threshold=0.01)
/* Calculates the cosine score between the two scans.
 * Arguments:
 *	first_scan - scan with mz:intensity.
 *	second_scan - scan with mz:intensity.
 *	threshold - mz tolerance to be considered the same peak.
 * Returns:
 *	score - the cosine score bewteen the two scans.
 */
{
	enforce(not_all_zeroes(first_scan.values),
			"The first scan is all 0, cosine is NaN");
	enforce(not_all_zeroes(second_scan.values),
			"The second scan is all 0, cosine is NaN");
	real[real][2] vectors = create_vectors(first_scan, second_scan, threshold);
	real[] first_scores = vectors[0].values;
	real[] second_scores = vectors[1].values;
	real score = calculate_cosine_score(first_scores, second_scores);
	return score;
}
unittest
{
	real[real] scan1 = [1: 2000, 2:3000, 3:1000];
	real[real] scan2 = [1: 2000, 2:3000, 3:1000];
	assert(approxEqual(find_cosine_score(scan1, scan2), 1));
	scan2 = [1: -2000, 2:-3000, 3:-1000];
	assert(approxEqual(find_cosine_score(scan1, scan2), -1));
	scan2 = [1: 0];
	assertThrown(find_cosine_score(scan1, scan2));
	assertThrown(find_cosine_score(scan2, scan1));
}
