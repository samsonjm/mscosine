/* Calculates cosine scores between MS scans.
 *
 * Author: Jonathan Samson
 * Date: 4 August 2020
 */
import std.algorithm;
import cosine;
import std.math;
import std.exception;
import std.string;
import std.conv;
import scans;
import std.stdio;
import std.getopt;
import mzxmlparser;

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

MS2Scan[] mgf_parser(string contents)
/* Parses the contents of a .mgf file into a list of MS2Scan objects.
 * Arguments:
 *	contents - the contents of a .mgf file.
 * Returns:
 *	scans - a list of MS2Scan objects populated by contents.
 */
{	
	MS2Scan[] scans;
	string[] lines = splitLines(contents);
	MS2Scan my_scan = new MS2Scan;
	foreach(string line; lines)
	{
		if (line == "BEGIN IONS")
			my_scan = new MS2Scan;
		else if (line == "END IONS")
			scans ~= my_scan;
		else if (line[0..1].isNumeric)
		{
			long space = indexOf(line, " ");
			my_scan.add_peak(line[0..space].to!real, 
					line[space+1..$].to!real);
		}
	}
	return scans;
}
unittest
{
	string contents = 
		"BEGIN IONS\nTITLE=12Conly_AA_MSMS_neg.1817.1817.1\nRTINSECONDS=647.13204\nPEPMASS=243.170868338055 253955.390625\nCHARGE=1-\n51.46782684 1460.6981201172\n75.82749939 1671.7169189453\n75.86730194 1605.3143310547\n100.1144104 1462.4990234375\n101.5387802 1490.517578125\n107.7608643 1808.1832275391\n118.443428 1619.8599853516\n130.0875244 37516.33203125\n146.9610138 1678.8117675781\n171.1526642 1760.8597412109\n199.1815948 35382.921875\n243.1713562 107272.828125\n244.1736908 8717.1875\nEND IONS\nBEGIN IONS\nTITLE=12Conly_AA_MSMS_neg.1825.1825.1\nRTINSECONDS=649.69068\nPEPMASS=243.170894922472 476078.830810500018\nCHARGE=1-\n94.02720642 1555.5870361328\n116.1556854 1499.4390869141\n129.1035614 3504.0708007813\n130.0875092 60504.5234375\n170.9055634 1636.6365966797\n199.1815796 54539.64453125\n200.1834564 1801.9443359375\n208.7902069 1524.8837890625\n213.2942657 2020.6203613281\n214.1232758 1974.2976074219\n243.171402 158421.265625\n244.1746063 11441.90625\nEND IONS";
	MS2Scan[] my_scans = mgf_parser(contents);
	assert(my_scans.length == 2);
	assert(my_scans[0].get_peak_intensity(101.5387802) == 
			1490.517578125);
	assert(my_scans[1].get_peak_intensity(101.5387802) != 
			1490.517578125);
}

void main(string[] args)
{
	string input_file;
	int scan_1_index;
	int scan_2_index;

	auto helpInformation = getopt(
			args,
			"input|i", "The input file in .mgl format", &input_file,
			"scan1|1", "The first scan to compare", &scan_1_index,
			"scan2|2", "The second scan to compare", &scan_2_index);
	if(helpInformation.helpWanted)
	{
		defaultGetoptFormatter(
				stdout.lockingTextWriter(),
				"Finds the cosine score between two scans.",
				helpInformation.options,
				"  %*s\t%*s%*s%s\n");
		return;
	}
	string file_contents = read_file(input_file);
	MS2Scan[] my_scans = mgf_parser(file_contents);
	real[real] peak_list_1 = my_scans[scan_1_index].get_peaks();
	real[real] peak_list_2 = my_scans[scan_2_index].get_peaks();

	real cosine_score = find_cosine_score(peak_list_1, peak_list_2);
	writeln(cosine_score);
}
