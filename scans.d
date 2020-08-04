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
	{
		return level;
	}

	real get_rt()
	{
		return retention_time;
	}

	real[real] get_peaks()
	{
		return peaks;
	}

	real get_peak_intensity(real my_peak)
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
	real parent_retention_time;

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

	void set_parent_rt(real rt)
	{
		parent_retention_time = rt;
	}

	Scan get_parent_scan()
	{
		return parent_scan;
	}

	real get_parent_peak()
	{
		return parent_peak;
	}

	real get_parent_rt()
	{
		return parent_retention_time;
	}
}
unittest
{
	Scan parent = new Scan;
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
	test.set_parent_rt(600.01);
	assert(test.get_parent_rt() == 600.01);
	test.set_parent_peak(244.1736908);
	assert(test.get_parent_peak() == 244.1736908);
}
