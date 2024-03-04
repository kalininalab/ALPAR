#ifndef SANS_TRANSLATOR_H
#define SANS_TRANSLATOR_H


#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include "util.h"

using namespace std;

/**
 * This namespace counts the occurrences where triplets have not been translated.
 */
namespace count {
    extern int getCount();
    extern void deleteCount();
}

/**
 * This class translates nucleotide sequences into amino acid sequences.
 */
class translator {

private:
    /**
     * This map contains the corresponding amino acid for each triplet.
     */
    static struct translationTable translationTable;

    /**
     * This set contains all allowed characters.
     */
    static unordered_set<char> bases;

    /**
     * This method checks if the given base triplet only contains allowed bases and contains 3 non-whitespace characters.
     * @param basicString a string which contains the next triplet
     */
    static bool checkUnit(string &basicString);

    /**
     * This method checks if the given char is a valid base character.
     * @param base the character
     */
    static bool isBase(char &base);
public:
    /**
     * This method initializes the translator.
     * @param codonfile the path to a valid file
     */
    static bool init(uint64_t& id);

    /**
     * This method returns the amino acid for the given triplet or - if any base is invalid - an empty string.
     * @param unit
     * @return amino acid
     */
    static string getTranslatedAminoAcid(string& unit);

    /**
     * This method takes a line of bases and translates them into amino acids by splitting them into triplets first.
     * @param line the line with bases
     * @return the translated line
     */
    static string translate(string& line);

protected:
    /**
     * This method initializes the data which is used to verify bases.
     */
    static void initAllowedBases();
};


#endif //SANS_TRANSLATOR_H