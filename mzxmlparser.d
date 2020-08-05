/* Tools to parse mzXML files into Scan[].
 * 
 * Author: Jonathan Samson
 * Date: 04-08-2020
 */
module mzxmldecoder;
import std.bitmanip;
import std.conv;
import std.base64;
import scans;
import std.stdio;
import std.math;
import std.exception;
import std.algorithm;

real[real] decode_mzxml_string(
		string encoded, 
		string compression="", 
		int precision=32)
{
	ubyte[] decoded = Base64.decode(encoded);
	real mz;
	real intensity;
	real[real] peak_list;
	int byte_size = precision / 8;

	if (compression=="zlib")
	{
		import std.zlib;
		decoded = cast(ubyte[]) uncompress(decoded);
	}
	for(int i = 1; i<=decoded.length/byte_size; ++i)
	{
		float readable;
		if (precision == 64)
		{
			ubyte[8] next_value = decoded[byte_size*(i-1)..
						     byte_size*i];
			readable = bigEndianToNative!double(next_value);
		}
		else // precision = 32
		{
			ubyte[4] next_value = decoded[byte_size*(i-1)..
						     byte_size*i];
			readable = bigEndianToNative!float(next_value);
		}

		if(i % 2 == 0)
		{
			intensity = readable.to!real;
			peak_list[mz] = intensity;
		}
		else
			mz = readable.to!real;
	}
	return peak_list;
}
unittest
{
	import std.algorithm;
	import std.math;
	real[real] answer = [
		51.4678:	1460.7,
		75.8275:	1671.72,
		75.8673:	1605.31,
		100.114:	1462.5,
		101.539:	1490.52,
		107.761:	1808.18,
		118.443:	1619.86,
		130.088:	37516.3,
		146.961:	1678.81,
		171.153:	1760.86,
		199.182:	35382.9,
		243.171:	107273,
		244.174:	8717.19
	];
	string line = "Qk3fDkS2lldCl6euRND28UKXvA9EyKoPQsg6lES2z/hCyxPbRLp" ~
		"QkELXhZBE4gXdQuzjCUTKe4VDAhZoRxKMVUMS9gVE0dn6QysnFUTcG4ND" ~
		"Ry59Rwo27ENzK95H0YRqQ3Qsd0YINMA=";
	real[real] function_test = decode_mzxml_string(line);
	assert(approxEqual(function_test.keys.sort, answer.keys.sort));
	assert(approxEqual(function_test.values.sort, answer.values.sort));

	line = "eJwBaACX/0JN3w5EtpZXQpenrkTQ9vFCl7wPRMiqD0LIOpREts/4QssT20" ~
		"S6UJBC14WQROIF3ULs4wlEynuFQwIWaEcSjFVDEvYFRNHZ+kMrJxVE3Bu" ~
		"DQ0cufUcKNuxDcyveR9GEakN0LHdGCDTAJ+wubA==";
	function_test = decode_mzxml_string(line, "zlib");
	assert(approxEqual(function_test.keys.sort, answer.keys.sort));
	assert(approxEqual(function_test.values.sort, answer.values.sort));

	line = "QEm74cAAAABAltLK4AAAAEBS9PXAAAAAQJoe3iAAAABAUveB4AAAAECZFU" ~
		"HgAAAAQFkHUoAAAABAltn/AAAAAEBZYntgAAAAQJdKEgAAAABAWvCyAAA" ~
		"AAECcQLugAAAAQF2cYSAAAABAmU9woAAAAEBgQs0AAAAAQOJRiqAAAABA" ~
		"Yl7AoAAAAECaOz9AAAAAQGVk4qAAAABAm4NwYAAAAEBo5c+gAAAAQOFG3" ~
		"YAAAABAbmV7wAAAAED6MI1AAAAAQG6FjuAAAABAwQaYAAAAAA==";
	function_test = decode_mzxml_string(line, "", 64);
	assert(approxEqual(function_test.keys.sort, answer.keys.sort));
	assert(approxEqual(function_test.values.sort, answer.values.sort));

	line = "eJxz8Nz98AADA4PDtEunHoDooC9fwfxZcvcUwPzvjWDxmaKOYDqSPagBrP" ~
		"7mfwYwP6k6AURP9xIC86M+bALTcxx2LwDRsXMSwebM9C8A8xOczoLlHwV" ~
		"2gflJcQfA9CxrewcQnZryCMyf3VwANjfj6Xkw/6HbXbC9eanVYPf9MugF" ~
		"q89r7QO76yDbDJC5AD9eO64=";
	function_test = decode_mzxml_string(line, "zlib", 64);
	assert(approxEqual(function_test.keys.sort, answer.keys.sort));
	assert(approxEqual(function_test.values.sort, answer.values.sort));
}

string read_file(string name_of_file)
/* Reads the file into a string.
 * Arguments:
 *      file_stream - The name of the file to read.
 *
 * Returns:
 *      file_contents - The contents of the file.
 */
{
        string file_contents = "";
        try
        {
                auto file = File(name_of_file, "r");
                string line;
                while ((line = file.readln()) !is null)
                {
                        file_contents ~= line;
                }
                file.close();
        }
        catch(ErrnoException e)
        {
                writeln("Invalid file name");
        }
        return file_contents;
}
unittest
{

}

Scan[] parse_mzxml(string contents)
/* Parses the contents of an .mzXML file into a list of Scan objects.
 * Arguments:
 *	contents - the contents of a .mzXML file.
 * Returns:
 *	scans - a list of Scan objects populated by contents.
 */
{
	Scan[] scans;
	return scans;
}
unittest
{
	string scans = read_file("example.mzXML");
	Scan[] parsed = parse_mzxml(scans);
	assert(approxEqual(parsed[0].get_rt(), 0.457129.to!real));
	assert(parsed[0].get_level() == 1);
	assert(approxEqual(parsed[2].get_rt(), 674.132.to!real));
	real[real] peaks = [
 		51.4678:        1460.7,
                75.8275:        1671.72,
                75.8673:        1605.31,
                100.114:        1462.5,
                101.539:        1490.52,
                107.761:        1808.18,
                118.443:        1619.86,
                130.088:        37516.3,
                146.961:        1678.81,
                171.153:        1760.86,
                199.182:        35382.9,
                243.171:        107273,
                244.174:        8717.19
	];
	assert(approxEqual(parsed[2].get_peaks().keys.sort, peaks.keys.sort));
	assert(approxEqual(parsed[2].get_peaks().values.sort,
			    peaks.values.sort));
	assert(parsed[2].get_level() == 2);
	assert(approxEqual(parsed[2].get_peak_intensity(171.153), 1760.86));
}
