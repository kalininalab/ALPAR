import os
import random
import string
import shutil


def generate_random_key():
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))

random_names = {}

def random_name_giver(strains_text_file, random_name_file_output):

    with open(strains_text_file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            strain = os.path.splitext(line.split("/")[-1].strip())[0]
            random_key = generate_random_key()
            # make sure it is unique
            while random_key in random_names.keys():
                random_key = generate_random_key()
            
            random_names[strain] = random_key

    with open(random_name_file_output, "w") as ofile:
        for key in random_names.keys():
            ofile.write(f"{key}\t{random_names[key]}\n")

    return random_names