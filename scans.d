/* Module for the scans class.
 *
 * Author: Jonathan Samson
 * Date: 04-08-2020
 */
module scans;

class Scan
{
	real retention_time;
	real[real] peaks;
	uint level;

	uint get_level()
	/* Gives the MS level of the scan.
	 * Returns:
	 *	this.level - The scan level.
	 */
	{
		return level;
	}

	real get_rt()
	/* Gives the retention time of the scan.
	 * Returns:
	 *	this.retention_time - The scan retention time.
	 */
	{
		return retention_time;
	}

	real[real] get_peaks()
	/* Gives the peaks of the scan.
	 * Returns:
	 *	this.peaks - The scan peaks as an associative array.
	 */
	{
		return peaks;
	}

	real get_peak_intensity(real my_peak)
	/* Gives the peak intensity of the set peak in thescan.
	 * Arguments:
	 *	my_peak - The peak of interest.
	 * Returns:
	 *	this.peaks[my_peak] - This scan's intensity of my_peak.
	 */
	{
		return peaks[my_peak];
	}

	void set_rt(real time)
	{
		retention_time = time;
	}

	void set_level(uint my_level)
	{
		level = my_level;
	}

	void set_peaks(real[real] my_peaks)
	{
		peaks = my_peaks;
	}

	void add_peak(real mz, real intensity)
	{
		peaks[mz] = intensity;
	}
}
unittest
{
	Scan test = new Scan;
	real[real] peaks = [
		51.46782684: 1460.6981201172,
		75.82749939: 1671.7169189453,
		75.86730194: 1605.3143310547,
		100.1144104: 1462.4990234375,
		101.5387802: 1490.517578125,
		107.7608643: 1808.1832275391,
		118.443428: 1619.8599853516,
		130.0875244: 37516.33203125,
		146.9610138: 1678.8117675781,
		171.1526642: 1760.8597412109,
		199.1815948: 35382.921875,
		243.1713562: 107272.828125,
		244.1736908: 8717.1875
	];
	test.set_level(1);
	assert(test.get_level() == 1);
	test.set_rt(100.110);
	assert(test.get_rt() == 100.110);
	test.set_peaks(peaks);
	assert(test.get_peaks() == peaks);
	test.add_peak(56.12356, 5235.12359);
	peaks[56.12356] = 5235.12359;
	assert(test.get_peak_intensity(56.12356) == 5235.12359);
	assert(test.get_peaks() == peaks);
}

class MS2Scan : Scan
{
	Scan parent_scan;
	real parent_peak;

	this()
	{
		level = 2;
	}

	void set_parent_scan(Scan my_parent)
	{
		parent_scan = my_parent;
	}

	void set_parent_peak(real peak) 
	{
		parent_peak = peak;
	}

	Scan get_parent_scan()
	/* Gives the parent scan of this MS2 scan.
	 * Returns:
	 *	this.parent_scan - The parent MS2 scan.
	 */
	{
		return parent_scan;
	}

	real get_parent_peak()
	/* Gives the peak from the parent scan that this MS2 refers to.
	 * Returns:
	 *	this.parent_peak - The relevant parent peak.
	 */
	{
		return parent_peak;
	}

	real get_parent_rt()
	/* Gives the retention time from the parent scan for this MS2.
	 * Returns:
	 *	this.parent_scan.get_rt() - The relevant parent rt.
	 */
	{
		return parent_scan.get_rt();
	}
}
unittest
{
	Scan parent = new Scan;
	parent.set_rt(600.100);
	Scan notparent = new Scan;
	MS2Scan test = new MS2Scan;
	real[real] peaks = [
		51.46782684: 1460.6981201172,
		75.82749939: 1671.7169189453,
		75.86730194: 1605.3143310547,
		100.1144104: 1462.4990234375,
		101.5387802: 1490.517578125,
		107.7608643: 1808.1832275391,
		118.443428: 1619.8599853516,
		130.0875244: 37516.33203125,
		146.9610138: 1678.8117675781,
		171.1526642: 1760.8597412109,
		199.1815948: 35382.921875,
		243.1713562: 107272.828125,
		244.1736908: 8717.1875
	];
	assert(test.get_level() == 2);
	test.set_rt(100.110);
	assert(test.get_rt() == 100.110);
	test.set_peaks(peaks);
	assert(test.get_peaks() == peaks);
	test.add_peak(56.12356, 5235.12359);
	peaks[56.12356] = 5235.12359;
	assert(test.get_peak_intensity(56.12356) == 5235.12359);
	assert(test.get_peaks() == peaks);
	test.set_parent_peak(244.1736908);
	assert(test.get_parent_peak() == 244.1736908);
	test.set_parent_scan(parent);
	assert(test.get_parent_scan() == parent);
	assert(test.get_parent_scan() != notparent);
	assert(test.get_parent_rt() == 600.100);
}
