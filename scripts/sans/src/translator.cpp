#include <sstream>
#include "translator.h"
#include "gc.h"

namespace count {
    int number = 0;

    int getCount() {
        return number;
    }

    void incrementCount() {
        number += 1;
    }

    void deleteCount() {
        number = 0;
    }
}

/**
 * codon and translationTable structs to save genetic code information
 */
struct codon {
    string triplet;
    string amino;
    bool isStart;
    bool isStop;

    codon() {
        triplet = "";
        amino = "";
        isStart = false;
        isStop = false;
    }

    bool isValid() const {
        return triplet.length() == 3 && amino.length() == 1;
    }

    bool operator==(const codon &other) const {
        return triplet == other.triplet;
    }

    size_t operator() (codon const& key) const {
        return hash<string>()(key.triplet);
    }
};
struct translationTable {

    static const char* bases[];
    unordered_set<string> names;
    uint64_t id;
    string ncbieaa;
    string sncbieaa;
    unordered_map<string, codon> codons;

    void createCodons(){
        if (!ncbieaa.empty() && !sncbieaa.empty()) {
            for (int i = 0; i <= 63; i++) {
                codon current = {};
                current.triplet += bases[0][i];
                current.triplet += bases[1][i];
                current.triplet += bases[2][i];
                current.amino += ncbieaa[i];
                current.isStop = sncbieaa[i] == '*';
                current.isStart = sncbieaa[i] == 'M';
                codons[current.triplet] = current;
            }
        }
    }

    codon operator[](const string& triplet) {
        auto it = codons.find(triplet);
        return it == codons.end() ? codon() :  it->second;
    }

};
struct translationTable translator::translationTable = {};
unordered_set<char> translator::bases;
const char* translationTable::bases[] = {"TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
                                         "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
                                         "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"};

bool translator::init(uint64_t& id) {
    string codonfile;

    for (int i = 0; i < config_gc_prt_len; i++) {
        codonfile += config_gc_prt[i];
    }

    std::stringstream ss(codonfile);
    std::string to;

    struct translationTable currentTranslationable = {};
    bool found = false;
    while (std::getline(ss,to,'\n') && !found) {
        util::trim(to);
        if (to.rfind("--", 0) != 0 && to.rfind("**", 0) != 0 && to.rfind("Genetic-code-table", 0) != 0 && !to.empty()) {

            util::replaceAll(to, "\"", "");
            util::replaceAll(to, ",", "");
            vector<string> parts = util::split(to, " ");

            if (to.rfind('{', 0) == 0) {
                currentTranslationable = {};
            } else if (parts.at(0) == "name") {
                currentTranslationable.names.insert(parts.at(1));
            } else if (parts.at(0) == "id") {
                currentTranslationable.id = stol(parts.at(1));
            } else if (parts.at(0) == "ncbieaa") {
                currentTranslationable.ncbieaa = parts.at(1);
            } else if (parts.at(0) == "sncbieaa") {
                currentTranslationable.sncbieaa = parts.at(1);
                currentTranslationable.createCodons();
            } else if (to.rfind('}', 0) == 0) {
                if (currentTranslationable.id == id) {
                    translator::translationTable = currentTranslationable;
                    found = true;
                }
            }
        }
    }
    if (!found) {
        cerr << "Error: Genetic code not found." << endl;
        exit(1);
    }
    return true;
}

void translator::initAllowedBases() {
    bases.insert('A');
    bases.insert('C');
    bases.insert('G');
    bases.insert('T');
}

string translator::getTranslatedAminoAcid(string &unit) {
    string translated;
    codon codon = translator::translationTable[unit];

    if (!codon.isValid()) {
        count::incrementCount();
        translated = "";
    } else {
        translated = codon.amino;
    }
    return  translated;
}

bool translator::checkUnit(string &basicString) {
    bool ok = true;

    // unit has to be 3 chars
    util::trim(basicString);
    ok = basicString.length() == 3;
    if (ok) {
        // only valid DNA bases are allowed
        ok = translator::isBase(basicString[0]);
        ok = ok && translator::isBase(basicString[1]);
        ok = ok && translator::isBase(basicString[2]);
    }

    return ok;
}

/**
 * This function checks if the character at the given position is allowed.
 * @param pos position in str
 * @param str the current part of the sequence
 * @return true if allowed, false otherwise
 */
bool translator::isBase(char &base) {
    unordered_set<char>::const_iterator got = translator::bases.find (base);
    return got != translator::bases.end();
}

string translator::translate(string &line) {
    string translated;
    for (auto pos = 0; pos+3 < line.length(); pos = pos +3) {
        string unit = line.substr(pos, 3);
        translated+= translator::getTranslatedAminoAcid(unit);
    }

    return translated;
}
