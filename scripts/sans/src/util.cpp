#include "util.h"


/**
 * This function compares the number of input genomes (n) to the compile parameter DmaxN (maxN).
 * Exits (code 3) if they disagree (too much) and a re-compilation is necessary (n>maxN) or recommended (n much smaller maxN).
 *
 * @param n number of input genomes
 * @return nothing
 */ 
void util::check_n(uint64_t& n, string& path, const uint64_t& max_N) {
	bool recompile = false;
	if (n>max_N) {
		cout << "Recompilation with DmaxN>=" << n << " necessary! (See " << path << "_autoN.)" << endl;
		recompile = true;
	} else if (max_N-n>=100) {
		cout << "Recompilation with DmaxN=" << n << " recommended. (See " << path << "_autoN.)" << endl;
		recompile = true;
	}
	if (recompile) {
		cout << "DmaxN\t" << max_N << endl;
		cout << "num\t" << n << endl;
		
		// open new makefile
		ofstream file(path+"_autoN");
		ostream newfilebuf(file.rdbuf());

		// open original make file
        ifstream infile(path);
        if (!infile.good()) {
            cerr << "Error: could not read makefile: "<< path << endl;
            exit(1);
        }
        
        // copy while replacing value for DmaxN and dependency on makefile
        regex ren("-DmaxN=[0-9]+");
        regex rem("makefile");
		string line,new_line;
        while (getline(infile, line)) {
			new_line=regex_replace(line, ren, "-DmaxN="+to_string(n));
			new_line=regex_replace(new_line, rem, "makefile_autoN");
			newfilebuf << new_line << endl;
		}
		
 		infile.close();
		file.close();
		
		// stop SANS
		exit(3);
	
	}
}




/**
 * This function calculates the arithmetic mean of two values.
 *
 * @param x first value
 * @param y second value
 * @return arithmetic mean
 */
double util::arithmetic_mean(uint32_t& x, uint32_t& y) {
    return x / 2.0 + y / 2.0;
}

/**
 * This function calculates the geometric mean of two values.
 *
 * @param x first value
 * @param y second value
 * @return geometric mean
 */
double util::geometric_mean(uint32_t& x, uint32_t& y) {
    return sqrt(x) * sqrt(y);
}

/**
 * This function calculates the geometric mean with pseudo counts.
 *
 * @param x first value
 * @param y second value
 * @return geometric mean
 */
double util::geometric_mean2(uint32_t& x, uint32_t& y) {
    return sqrt(x+1) * sqrt(y+1);
}

/**
 * This function displays a duration in a human readable format.
 *
 * @param time duration
 * @return formatted string
 */
string util::format_time(chrono::high_resolution_clock::duration time) {
    double value = chrono::duration_cast<chrono::milliseconds>(time).count();
    string unit = "ms";
    if (value >= 1000) {
        value /= 1000;
        unit = "sec";
        if (value >= 200) {
            value /= 60;
            unit = "min";
            if (value >= 200) {
                value /= 60;
                unit = "h";
            }
        }
    }
    string number = to_string(value);
    return number.substr(0, number.find('.')+2) + ' ' + unit;
}

/**
* Char to bit and bit to char transcoding functions  
*/
/**
 * This function encodes a single character to two bits.
 *
 * @param c character
 * @return bit sequence
 */
uint64_t util::char_to_bits(const char& c) {
    switch (c) {
        case 'A':
            return 0b00u;
        case 'C':
            return 0b01u;
        case 'G':
            return 0b10u;
        case 'T':
            return 0b11u;
        default:
            cerr << "Error: Invalid character " << c << "." << endl;
            return -1;
    }
}

/**
 * This function decodes two bits to a single character.
 *
 * @param b bit sequence
 * @return character
 */
char util::bits_to_char(uint64_t& b) {
    switch (b) {
        case 0b00u:
            return 'A';
        case 0b01u:
            return 'C';
        case 0b10u:
            return 'G';
        case 0b11u:
            return 'T';
        default:
            cerr << "Error: Invalid bit sequence " << b << "." << endl;
            return -1;
    }
}


/**
 * This function encodes a single character to five bits.
 *
 * @param c character
 * @return bit sequence
 */
uint64_t util::amino_char_to_bits(char& c) {

    switch (c) {
        case 'A':
            return 0b00000u;
        case 'B':
            return 0b00001u;
        case 'C':
            return 0b00010u;
        case 'D':
            return 0b00011u;
        case 'E':
            return 0b00100u;
        case 'F':
            return 0b00101u;
        case 'G':
            return 0b00110u;
        case 'H':
            return 0b00111u;
        case 'I':
            return 0b01000u;
        case 'J':
            return 0b01001u;
        case 'K':
            return 0b01010u;
        case 'L':
            return 0b01011u;
        case 'M':
            return 0b01100u;
        case 'N':
            return 0b01101u;
        case 'O':
            return 0b01110u;
        case 'P':
            return 0b01111u;
        case 'Q':
            return 0b10000u;
        case 'R':
            return 0b10001u;
        case 'S':
            return 0b10010u;
        case 'T':
            return 0b10011u;
        case 'U':
            return 0b10100u;
        case 'V':
            return 0b10101u;
        case 'W':
            return 0b10110u;
        case 'X':
            return 0b10111u;
        case 'Y':
            return 0b11000u;
        case 'Z':
            return 0b11001u;
        case '*':
            return 0b11010u;
        default:
            cerr << "Error: Invalid character " << c << "." << endl;
            exit(1);
    }
}

/**
 * This function decodes five bits to a single character.
 *
 * @param b bit sequence
 * @return character
 */
char util::amino_bits_to_char(uint64_t& b) {
    switch (b) {
        case 0b00000u:
            return 'A';
        case 0b00001u:
            return 'B';
        case 0b00010u:
            return 'C';
        case 0b00011u:
            return 'D';
        case 0b00100u:
            return 'E';
        case 0b00101u:
            return 'F';
        case 0b00110u:
            return 'G';
        case 0b00111u:
            return 'H';
        case 0b01000u:
            return 'I';
        case 0b01001u:
            return 'J';
        case 0b01010u:
            return 'K';
        case 0b01011u:
            return 'L';
        case 0b01100u:
            return 'M';
        case 0b01101u:
            return 'N';
        case 0b01110u:
            return 'O';
        case 0b01111u:
            return 'P';
        case 0b10000u:
            return 'Q';
        case 0b10001u:
            return 'R';
        case 0b10010u:
            return 'S';
        case 0b10011u:
            return 'T';
        case 0b10100u:
            return 'U';
        case 0b10101u:
            return 'V';
        case 0b10110u:
            return 'W';
        case 0b10111u:
            return 'X';
        case 0b11000u:
            return 'Y';
        case 0b11001u:
            return 'Z';
        case 0b11010u:
            return '*';
        default:
            cerr << "Error: Invalid bit sequence " << b << "." << endl;
            exit(1);
    }
}

void util::ltrim(std::string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

void util::rtrim(std::string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

void util::trim(std::string &s) {
    util::ltrim(s);
    util::rtrim(s);
}

void util::replaceAll(string &s, string replaceString, string newString){
    if (!replaceString.empty()) {
        size_t start_pos = 0;
        while ((start_pos = s.find(replaceString, start_pos)) != std::string::npos) {
            s.replace(start_pos, replaceString.length(), newString);
            start_pos += newString.length();
        }
    }
}

vector<string> util::split(string &s, string delimiter) {
    vector<string> vector;
    auto start = 0U;
    auto end = s.find(delimiter);
    string basicString;
    while (end != std::string::npos) {
        basicString = s.substr(start, end - start);
        if (!basicString.empty()) {
            vector.push_back(basicString);
        }
        start = end + delimiter.length();
        end = s.find(delimiter, start);
    }

    basicString = s.substr(start, end);
    if (!basicString.empty()) {
        vector.push_back(basicString);
    }

    return vector;
}

bool util::is_number(string& s) {
    return !s.empty() && find_if(s.begin(), s.end(), [](unsigned char c) { return !isdigit(c); }) == s.end();
}



bool util::path_exist(const string &filename) {
	string folder="";
	uint64_t found=filename.find_last_of("/\\");
	if (found!=string::npos)
    { 
        folder=filename.substr(0,found+1);
		struct stat buffer;
		return (stat (folder.c_str(), &buffer) == 0);
    } else {
		return true;
	}
}


