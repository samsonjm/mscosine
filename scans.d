/* Module for the scans class.
 *
 * Author: Jonathan Samson
 * Date: 04-08-2020
 */
module scans;

class Scan
{
	ureal retention_time;
	ureal[ureal] peaks;
	uint level;
}

class MS2Scan : Scan
{
	level = 2;
	*Scan parent_scan;
	ureal parent_peak;
	ureal parent_retention_time;
}
