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
}

class MS2Scan : Scan
{
	Scan *parent_scan;
	real parent_peak;
	real parent_retention_time;

	this()
	{
		level = 2;
	}
}
