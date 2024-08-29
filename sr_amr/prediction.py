import os
import pickle
import csv
import numpy as np

def process_data_for_prediction(input_path_or_file):

    if os.path.isdir(input_path_or_file):
        input_dir = input_path_or_file
        input_files = os.listdir(input_dir)
    
    else:
        input_file = input_path_or_file
        with open(input_file, 'r') as infile:
            input_files = infile.readlines()
    
    return input_files

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