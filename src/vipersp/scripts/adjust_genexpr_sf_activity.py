import argparse
import torch
import pandas as pd
import numpy as np
from vipersp.model import EWlayer, FClayer

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}

##### FUNCTIONS #####
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--activity_file", type=str)
    parser.add_argument("--weights_files", type=str)
    parser.add_argument("--input_regulators_file", type=str)
    parser.add_argument("--output_regulators_file", type=str)
    parser.add_argument("--model_type", type=str)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args

def main():
    args = parse_args()
    activity_file = args.activity_file
    weights_files = args.weights_files
    input_regulators_file = args.input_regulators_file
    output_regulators_file = args.output_regulators_file
    model_type = args.model_type
    output_file = args.output_file

    # load
    activity = pd.read_table(activity_file, index_col=0)
    weights = [torch.load(f) for f in weights_files.split(",")]
    input_regulators = list(pd.read_table(input_regulators_file, header=None)[0])
    output_regulators = list(pd.read_table(output_regulators_file, header=None)[0])

    # prep data
    activity = pd.merge(
        pd.DataFrame(index=input_regulators),
        activity,
        how="left", left_index=True, right_index=True
    )
    X = torch.tensor(activity.fillna(0).T.values, dtype=torch.float32)

    if model_type=="fclayer":
        model = FClayer(input_size=len(input_regulators), output_size=len(output_regulators))

    elif model_type=="ewlayer":
        model = EWlayer(input_size=len(input_regulators))

    activity_preds = []
    for k in range(len(weights)):
        model.load_state_dict(weights[k])        
        # make predictions
        model.eval()
        with torch.no_grad():
            Y_hat = model(X)

        # prep outputs
        activity_pred = pd.DataFrame(Y_hat.detach().numpy().T, index=output_regulators, columns=activity.columns)
        activity_pred.index.name = "regulator"
        
        activity_preds.append(activity_pred)
        
    # average predictions
    activity_stack = np.stack([df.values for df in activity_preds])
    activity_avg = np.mean(activity_stack, axis=0)
    activity_preds = pd.DataFrame(
        activity_avg, index=activity_preds[0].index, columns=activity_preds[0].columns
    )
    
    # save
    activity_preds.reset_index().to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")