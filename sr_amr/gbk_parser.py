from Bio import SeqIO
import os

class Records:
    """
    Subclass for store the record information
    """
    def __init__(self):
        self.gene = None
        self.start = None
        self.end = None
        self.translation = None
        self.protein_id = None
        self.organism = None
        self.locus_tag = None
        self.EC_Number = None
        self.interface = None
        self.codon_start = None
        self.transl_table = None
        self.product = None
        self.db_xref = None
        self.record_id = None
        self.feature_type = None
        self.mol_type = None
        self.protein_gpa_name = None

def gbk_parser(input_strain_name, input_gbk):
    """
    Gets gbk file as input, parses it and returns list of Records class

    Parameters
    ----------
    input_gbk : str
        path to gbk file

    Return
    ----------
    all_the_records : list
        list of Records class

    """

    recs = [rec for rec in SeqIO.parse(input_gbk, "genbank")]

    all_the_records = []

    for rec in recs:

        feats = [feat for feat in rec.features]

        for feat in feats:

            temp = Records()

            if 'gene' in feat.qualifiers.keys():
                temp.gene = feat.qualifiers['gene'][0]

            temp.start = int(feat.location.start)
            
            temp.end = int(feat.location.end)

            temp.record_id = rec.id
            
            temp.feature_type = feat.type

            if 'feature_type' in feat.qualifiers.keys():
                temp.translation = feat.qualifiers['translation'][0]

            if 'protein_id' in feat.qualifiers.keys():
                temp.protein_id = feat.qualifiers['protein_id'][0]
            
            if 'organism' in feat.qualifiers.keys():
                temp.organism = feat.qualifiers['organism'][0]
            
            if 'locus_tag' in feat.qualifiers.keys():
                temp.locus_tag = feat.qualifiers['locus_tag'][0]

            if 'EC_Number' in feat.qualifiers.keys():
                temp.EC_Number = feat.qualifiers['EC_Number'][0]

            if 'interface' in feat.qualifiers.keys():
                temp.interface = feat.qualifiers['interface'][0]

            if 'codon_start' in feat.qualifiers.keys():
                temp.codon_start = feat.qualifiers['codon_start'][0]

            if 'transl_table' in feat.qualifiers.keys():
                temp.transl_table = feat.qualifiers['transl_table'][0]

            if 'product' in feat.qualifiers.keys():
                temp.product = feat.qualifiers['product'][0]

            if 'db_xref' in feat.qualifiers.keys():
                temp.db_xref = feat.qualifiers['db_xref'][0]

            if 'mol_type' in feat.qualifiers.keys():
                temp.mol_type = feat.qualifiers['mol_type'][0]
            
            if temp.locus_tag and temp.product:
                temp.protein_gpa_name = f"{input_strain_name};{temp.locus_tag};{temp.product}"

            all_the_records.append(temp)

    return all_the_records

# Usage example;
#test_results = gbk_parser(f"{os.path.dirname(os.path.realpath(__file__))}/PROKKA_05262023.gbk")