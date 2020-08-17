import std.stdio;
import scans;
import cosine;
import mscosine;
import mzxmlparser;
import std.getopt;
import std.regex;

void main(string[] args)
{
       string input_file;
        int scan_1_index;
        int scan_2_index;

        auto helpInformation = getopt(
                        args,
                        "input|i", "The input file in .mgl or .mzxml format",
                        &input_file,
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
        auto file_extension = ctRegex!(`\.(\w*)$`);
        MSXScan[] my_scans;
        switch (input_file.matchFirst(file_extension)[1])
        {
                default:
                {
                        throw new Exception("Invalid input file extension.");
                }
                case "mgl":
                {
                       my_scans = mgf_parser(file_contents);
                        break;
                }
                case "mzXML":
                {
                        my_scans = parse_mzxml(file_contents);
                        break;
                }
        }
        real[real] peak_list_1 = my_scans[scan_1_index-1].peaks;
        real[real] peak_list_2 = my_scans[scan_2_index-1].peaks;
        real cosine_score = find_cosine_score(peak_list_1, peak_list_2);
        writeln(cosine_score);
}
