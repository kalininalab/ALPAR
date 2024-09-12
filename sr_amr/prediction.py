import os
import pickle
import csv
import numpy as np
from sr_amr.ml import prps_ml_preprecessor

def process_data_for_prediction(input_path_or_file):

    if os.path.isdir(input_path_or_file):
        input_dir = input_path_or_file
        input_files = os.listdir(input_dir)
    
    else:
        input_file = input_path_or_file
        with open(input_file, 'r') as infile:
            input_files = infile.readlines()
    
    return input_files


def equalize_columns(binary_table1, binary_table2, output_file):
    # Read binary_table1 into a dictionary of dictionaries
    with open(binary_table1, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        headers1 = next(reader)
        data1 = {rows[0]: {headers1[i]: rows[i] for i in range(1, len(headers1))} for rows in reader}

    # Read binary_table2 into a dictionary of dictionaries
    with open(binary_table2, 'r') as infile2:
        reader = csv.reader(infile2, delimiter='\t')
        headers2 = next(reader)
        data2 = {rows[0]: {headers2[i]: rows[i] for i in range(1, len(headers2))} for rows in reader}

    # Find columns that are in binary_table1 but not in binary_table2
    missing_columns = [col for col in headers1 if col not in headers2]

    # Find columns that are in binary_table2 but not in binary_table1
    extra_columns = [col for col in headers2 if col not in headers1]

    # Add missing columns to binary_table2 with default value '0'
    for key in data2:
        for col in missing_columns:
            data2[key][col] = '0'

    # Remove extra columns from binary_table2
    for key in data2:
        for col in extra_columns:
            data2[key].pop(col)

    # Reorder columns in data2 to match the order in headers1
    headers2 = headers1
    reordered_data2 = {}
    for key in data2:
        reordered_data2[key] = {col: data2[key].get(col, '0') for col in headers1 if col not in ['Strain', '', 'strain', 'Strains', 'strains', ' ']}

    # Write the modified binary_table2 to the output file
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers1)
        for key, value in reordered_data2.items():
            writer.writerow([key] + [value[col] for col in headers1 if col not in ['Strain', '', 'strain', 'Strains', 'strains', ' ']])


def predict(trained_model, prediction_df, output_dir):

    prediction_dict = {}

    # Load the trained model
    with open(trained_model, 'rb') as trained_model_file:
        model = pickle.load(trained_model_file)

    # Read the CSV file into a list of lists
    with open(prediction_df, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        data = list(reader)

    # Extract the index (first column) and the genotype data
    indices = [row[0] for row in data]
    genotype_array = np.array([row[1:] for row in data], dtype=float)

    # Use the loaded model to make predictions
    predictions = model.predict(genotype_array)

    # Create the prediction dictionary
    for idx, prediction in zip(indices, predictions):
        prediction_dict[idx] = prediction

    # Write the predictions to a file
    with open(os.path.join(output_dir, "predictions.csv"), "w") as outfile:
        for key in prediction_dict.keys():
            outfile.write(f"{key}\t{prediction_dict[key]}\n")

# Example usage
# predict('path_to_trained_model.pkl', 'path_to_prediction_df.tsv', 'output_directory')