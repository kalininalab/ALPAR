import os

def fia_file_annotation(fia_file, annotation_file):
    with open (fia_file, "r") as ifile:
        lines = ifile.readlines()
    
    parent_dir = os.path.dirname(fia_file)
    filename = os.path.basename(fia_file)

    with open(annotation_file, "r") as annotation_file_handle:
        annotation_lines = annotation_file_handle.readlines()

    annotation_dict = {}

    for annotation_line in annotation_lines:
        splitted = annotation_line.split("\t")
        annotation_dict[splitted[0].strip()] = [splitted[1].strip(), splitted[2].strip(), splitted[3].strip()]

    with open(os.path.join(parent_dir, f"Annotated_{filename}"), "w") as ofile:
        ofile.write("Mutation;ImportanceValue;Effect;Gene;Product\n")
        for line in lines:
            splitted = line.split(";")
            mutation = splitted[0].strip()
            importance = splitted[1].strip()
            effect = annotation_dict[mutation][0].strip()
            gene = annotation_dict[mutation][1].strip()
            product = annotation_dict[mutation][2].strip()
            ofile.write(f"{mutation};{importance};{effect};{gene};{product}\n")