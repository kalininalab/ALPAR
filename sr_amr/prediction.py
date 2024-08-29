import os
import pickle
import pandas as pd
import numpy as np

def predict(trained_model, prediction_df, output_dir):

    prediction_dict = {}

    with open(trained_model, 'rb') as trained_model_file:
        model = pickle.load(trained_model_file)


    genotype_df = pd.read_csv(prediction_df, sep="\t", index_col=0, header=0)
    genotype_array = genotype_df.to_numpy()

    #print(genotype_df)
    #print(genotype_array)

    test_data = np.array(genotype_array)

    # Assuming test_data is your test data
    # Use the loaded model to make predictions
    predictions = model.predict(genotype_array)
    #predictions_proba = model.predict_proba(genotype_array)

    # cnt = 0
    # for prediction in predictions_proba:
    #     prediction_dict[genotype_df.index[cnt]] = prediction
    #     cnt += 1

    cnt = 0
    for prediction in predictions:
        prediction_dict[genotype_df.index[cnt]] = prediction
        cnt += 1

    with open(os.path.join(output_dir, "predictions.csv"), "w") as outfile:
        for key in prediction_dict.keys():
            outfile.write(f"{key}\t{prediction_dict[key]}\n")

predict()