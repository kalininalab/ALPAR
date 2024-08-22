import os
from Bio import SeqIO
import textwrap
import warnings

warnings.filterwarnings("ignore")

def annotation_file_creator(reference_genome, output_folder):

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


    recs = [rec for rec in SeqIO.parse(reference_genome, "genbank")]

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

            all_the_records.append(temp)

    all_the_records = sorted(all_the_records, key=lambda record: record.start)

    with open(os.path.join(output_folder, "annotation_file.tsv"), "w") as ofile:
        ofile.write("record_id\tgene\tstart\tend\ttranslation\tprotein_id\torganism\tlocus_tag\tEC_Number\tinterface\tcodon_start\ttransl_table\tproduct\tdb_xref\tfeature_type\tmol_type\n")
        for record in all_the_records:
            ofile.write(f"{record.record_id}\t{record.gene}\t{record.start}\t{record.end}\t{record.translation}\t{record.protein_id}\t{record.organism}\t{record.locus_tag}\t{record.EC_Number}\t{record.interface}\t{record.codon_start}\t{record.transl_table}\t{record.product}\t{record.db_xref}\t{record.feature_type}\t{record.mol_type}\n")

    return all_the_records

# mutations_annotation function should transfer to the snippy output check annotation_file_from_snippy function
def mutations_annotation_adder(annotation_file ,output_folder, all_the_records):

    mutations_position_dict = {}

    with open(annotation_file, "r") as infile:
        annotation_lines = infile.readlines()
        for header_line in annotation_lines:
            header_line_splitted = header_line.split("\t")
            mutation = header_line_splitted[0]
            position = mutation.split(",")[0].strip()[1:-1].strip()
            mutations_position_dict[mutation] = position
        
    # Sort mutations_position_dict by values
    mutations_position_dict = {k: v for k, v in sorted(mutations_position_dict.items(), key=lambda item: item[1])}

    with open(os.path.join(output_folder, "mutations_annotations.tsv"), "w") as ofile:
        for mutation in mutations_position_dict.keys():
            for record in all_the_records:
                position = mutations_position_dict[mutation]
                if int(record.start) <= int(position) <= int(record.end):
                    ofile.write(f"{mutation}\t{record.start}\t{record.end}\t{record.gene}\t{record.product}\t{record.translation}\t{record.protein_id}\t{record.db_xref}\n")
                    break

def annotation_function(annotation_file, output_folder, reference_genome):

    all_the_records = annotation_file_creator(reference_genome, output_folder)

    mutations_annotation_adder(annotation_file, output_folder, all_the_records)


def structman_input_creator(args):
        
    is_folder = False

    if os.path.isdir(args.input):
        is_folder = True

    amino_acid_three_letter_to_one_letter_dict = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
                                                  'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
                                                  'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
                                                  'Ala': 'A', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}

    amino_acid_one_letter_to_three_letter_dict = {'C': 'Cys', 'D': 'Asp', 'S': 'Ser', 'Q': 'Gln', 'K': 'Lys', 'I': 'Ile', 'P': 'Pro', 'T': 'Thr',
                                                  'F': 'Phe', 'N': 'Asn', 'G': 'Gly', 'H': 'His', 'L': 'Leu', 'R': 'Arg', 'W': 'Trp', 'A': 'Ala', 'V': 'Val', 'E': 'Glu', 'Y': 'Tyr', 'M': 'Met'}

    if is_folder:
        # paths will be corrected after getting first results
        translate_dict = {}
        annotation_function(args.annotation, os.path.join(
            args.temp, "structman"), args.reference_genome)
        created_annotation_file = os.path.join(
            args.temp, "structman", "mutations_annotations.tsv")
        with open(created_annotation_file, "r") as annot_file:
            annot_file_lines = annot_file.readlines()
            for annot_file_line in annot_file_lines:
                splitted = annot_file_line.split("\t")
                translate_dict[splitted[0].strip()] = splitted[5].strip()

        positions_dict = {}
        gene_dict = {}
        with open(args.annotation, "r") as positions_file:
            position_lines = positions_file.readlines()
            for position_line in position_lines[1:]:
                splitted = position_line.split("\t")
                mutation = splitted[0].strip()
                mutation_type = mutation.split(",")[2].strip()[1:-1]
                # so far only support snps, complexes will be add in the future
                if mutation_type == "snp":
                    effect = splitted[1].strip()
                    ref_aa = effect.split(" ")[-1][2:5]
                    alt_aa = effect.split(" ")[-1][-3:]
                    pos_aa = effect.split(" ")[-1][5:-3]
                    positions_dict[mutation] = [ref_aa, int(pos_aa), alt_aa]
                    if splitted[2].strip() == "":
                        gene_dict[mutation] = "unknown"
                    else:
                        gene_dict[mutation] = splitted[2].strip()


        if os.path.exists(os.join(args.input, "ml")):
            ml_outputs = os.listdir(os.join(args.input, "ml"))
            with open(os.path.join(args.output, "structman", "structman_input.smlf"), "w") as structman_input:
                for ml_output in ml_outputs:
                    if ml_output.startswith("Annotated_"):
                        annotated_ml_output = os.join(args.input, "ml", ml_output)
                        with open(annotated_ml_output, "r") as infile:
                            lines = infile.readlines()
                            if len(lines) > 0:
                                for line in lines:
                                    if line != "":
                                        splitted = line.split(";")
                                        mutation = splitted[0].strip()
                                        if (mutation in translate_dict) and (mutation in positions_dict) and (mutation in gene_dict):
                                            header_line = f">{args.name}:{gene_dict[mutation]} {amino_acid_three_letter_to_one_letter_dict[positions_dict[mutation][0]]}{positions_dict[mutation][1]}{amino_acid_three_letter_to_one_letter_dict[positions_dict[mutation][2]]}"
                                            translation_line = f"{translate_dict[mutation]}"
                                            translation_line_wrapped = "\n".join(textwrap.wrap(translation_line, width=70))
                                            structman_input.write(header_line + "\n" + translation_line_wrapped + "\n")