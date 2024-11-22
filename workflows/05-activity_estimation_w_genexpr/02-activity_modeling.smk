import os
import copy

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","activity_estimation_w_genexpr")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

DATASETS = {
    "CCLE": {
        "genexpr": os.path.join(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz"), 
        "EX": os.path.join(PREP_DIR,"event_psi","CCLE-EX.tsv.gz")
    }
}

OMIC_PERT_DICT = {
    "EX": "EX",
    "bulkgenexpr": "genexpr",
    "scgenexpr": "genexpr",
    "bulkscgenexpr": "genexpr"
}

REGULON_DIRS = {
    "bulkgenexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-bulkgenexpr"),
    "scgenexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-scgenexpr"),
    "bulkscgenexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-bulkscgenexpr"),
    "EX": os.path.join(ROOT,"results","new_empirical_network","files","experimentally_derived_regulons_pruned_w_viper_networks-EX")
}

METHODS_ACTIVITY = ["viper","correlation_pearson","correlation_spearman","gsea"]

OMIC_TYPES = ["bulkgenexpr","scgenexpr","bulkscgenexpr","EX"]
OMIC_GENEXPR_REGULONS = ["bulkgenexpr","scgenexpr","bulkscgenexpr"]
MODEL_ARCHS = ["fclayer","ewlayer"]

##### RULES #####
rule all:
    input:
        # signature CCLE
        expand(os.path.join(RESULTS_DIR,"files","signatures","CCLE-{omic_signature}.tsv.gz"), omic_signature=["EX","genexpr"]),

        # estimate splicing factor activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","CCLE-{omic_regulon}.tsv.gz"), omic_regulon=OMIC_TYPES),
        
        # train models
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","losses.tsv.gz"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_ARCHS),
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights.pth"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_ARCHS),
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","common_regulators.tsv.gz"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_ARCHS),
        os.path.join(RESULTS_DIR,"files","model_sf_activity","losses-merged.tsv.gz")
        
        # # estimate SF activities with trained models
        # expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_{model_type}_and_{omic_regulon}.tsv.gz"), dataset=PERT_SPLICING_FILES.keys(), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_ARCHS),

        # make figures
        # os.path.join(RESULTS_DIR,"figures","activity_modeling")
        
        
rule compute_signature_within:
    input:
        data = lambda wildcards: DATASETS["CCLE"][wildcards.omic_signature]
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","CCLE-{omic_signature}.tsv.gz")
    run:
        import pandas as pd
        
        data = pd.read_table(input.data, index_col=0)
        
        # subtract median
        signature = data
        #signature = signature - signature["ACH-001086"].values.reshape(-1,1)
        signature = signature - signature.median(axis=1).values.reshape(-1,1)
        
        # save
        signature.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")

rule compute_protein_activity:
    input:
        signature = lambda wildcards: os.path.join(RESULTS_DIR,"files","signatures","CCLE-{omic_signature}.tsv.gz").format(omic_signature=OMIC_PERT_DICT[wildcards.omic_regulon]),
        regulons_path = lambda wildcards: REGULON_DIRS[wildcards.omic_regulon]
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","CCLE-{omic_regulon}.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
        """
        
rule train_model_to_adjust_sf_activity:
    input:
        features = os.path.join(RESULTS_DIR,"files","protein_activity","CCLE-{omic_regulon}.tsv.gz"),
        labels = os.path.join(RESULTS_DIR,"files","protein_activity","CCLE-EX.tsv.gz")
    output:
        losses = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","losses.tsv.gz"),
        weights = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights.pth"),
        common_regulators = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","common_regulators.tsv.gz")
    params:
        model_type = "{model_type}",
        batch_size = 128,
        max_epochs = 20,
        learning_rate = 0.01,
        train_split = 0.85
    threads: 4
    run:
        import numpy as np
        import pandas as pd
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset
        import lightning as L
        from vipersp.model import EWlayer, FClayer, LitWrapper
        L.seed_everything(1234, workers=True)
        
        # load data
        features = pd.read_table(input.features, index_col=0)
        labels = pd.read_table(input.labels, index_col=0)

        # unpack params
        model_type = params.model_type
        batch_size = params.batch_size
        max_epochs = params.max_epochs
        learning_rate = params.learning_rate
        train_split = params.train_split
        num_workers = threads
        
        # split
        avail_samples = list(set(features.columns).intersection(labels.columns))
        train_samples = np.random.choice(avail_samples, int(len(avail_samples)*train_split))
        val_samples = list(set(avail_samples) - set(train_samples))
        
        common_regulators = list(set(features.index).intersection(labels.index))
        
        train_features = features.loc[common_regulators, train_samples].copy()
        train_labels = labels.loc[common_regulators, train_samples].copy()
        val_features = features.loc[common_regulators, val_samples].copy()
        val_labels = labels.loc[common_regulators, val_samples].copy()
        
        # prep dataloaders
        X_train = torch.tensor(train_features.fillna(0).T.values, dtype=torch.float32)
        Y_train = torch.tensor(train_labels.fillna(0).T.values, dtype=torch.float32)
        X_val = torch.tensor(val_features.fillna(0).T.values, dtype=torch.float32)
        Y_val = torch.tensor(val_labels.fillna(0).T.values, dtype=torch.float32)
        
        ds_train = TensorDataset(X_train, Y_train)
        ds_val = TensorDataset(X_val, Y_val)
        
        train_loader = DataLoader(ds_train, batch_size=batch_size, shuffle=True, num_workers=num_workers)
        val_loader = DataLoader(ds_val, batch_size=batch_size, shuffle=False, num_workers=num_workers)
        
        # train model
        criterion = nn.SmoothL1Loss()
        if model_type=="fclayer":
            model = FClayer(input_size=len(common_regulators), output_size=len(common_regulators))
        
        elif model_type=="ewlayer":
            model = EWlayer(input_size=len(common_regulators))        
            
        litwrapper = LitWrapper(model, criterion, learning_rate)
        trainer = L.Trainer(
            max_epochs = max_epochs,
            accelerator = "cpu",
            enable_checkpointing = False, # do not save checkpoints
            logger = False # do not save logs
        )

        # model iterations are being saved as we train
        trainer.fit(litwrapper, train_loader, val_loader)
        
        # prepare outputs
        losses = pd.DataFrame(litwrapper.losses)
        losses = losses.drop(columns="batch").groupby(["epoch"]).mean().reset_index()
        losses = losses.melt(id_vars="epoch", value_vars=["val_loss","train_loss","val_pearson","train_pearson"], var_name="loss_type", value_name="loss")
        
        # save
        ## loss
        losses.to_csv(output.losses, **SAVE_PARAMS)
        ## model weights        
        torch.save(model.state_dict(), output.weights)
        ## common regulators
        pd.DataFrame(common_regulators).to_csv(output.common_regulators, header=None, **SAVE_PARAMS)
        
        print("Done!")
        

rule combine_losses:
    input:
        losses = [os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","losses.tsv.gz").format(omic_regulon=o, model_type=m) for o in OMIC_GENEXPR_REGULONS for m in MODEL_ARCHS]
    output:
        loss = os.path.join(RESULTS_DIR,"files","model_sf_activity","losses-merged.tsv.gz")
    run:
        import os
        import pandas as pd
        
        losses = []
        for f in input.losses:
            loss = pd.read_table(f)
            loss["model_type"] = os.path.basename(os.path.dirname(f))
            loss["omic_regulon"] = os.path.basename(os.path.dirname(os.path.dirname(f)))
            losses.append(loss)
            
        losses = pd.concat(losses)
        
        losses.to_csv(output.loss, **SAVE_PARAMS)
        
        print("Done!")
        
        
# rule predict_sf_activity_from_model:
#     input:
#         activity = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz"),
#         weights = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights.pth"),
#         common_regulators = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","common_regulators.tsv.gz")
#     output:
#         activity_pred = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_{model_type}_and_{omic_regulon}.tsv.gz")
#     params:
#         model_type = "{model_type}"        
#     run:
#         import torch
#         import pandas as pd
#         from vipersp.model import EWlayer, FClayer
        
#         # load
#         activity = pd.read_table(input.activity, index_col=0)
#         weights = torch.load(input.weights)
#         common_regulators = list(pd.read_table(input.common_regulators, header=None)[0])
#         model_type = params.model_type
        
#         # prep data
#         activity = pd.merge(
#             pd.DataFrame(index=common_regulators),
#             activity,
#             how="left", left_index=True, right_index=True
#         )
#         X = torch.tensor(activity.fillna(0).T.values, dtype=torch.float32)
        
#         if model_type=="fclayer":
#             model = FClayer(input_size=len(common_regulators), output_size=len(common_regulators))
        
#         elif model_type=="ewlayer":
#             model = EWlayer(input_size=len(common_regulators))        
            
#         model.load_state_dict(weights)        
        
#         # make predictions
#         model.eval()
#         with torch.no_grad():
#             Y_hat = model(X)
        
#         # prep outputs
#         activity_pred = pd.DataFrame(Y_hat.detach().numpy().T, index=common_regulators, columns=activity.columns)
#         activity_pred.index.name = "regulator"
        
#         # save
#         activity_pred.reset_index().to_csv(output.activity_pred, **SAVE_PARAMS)
        
#         print("Done!")
        
        
# rule figures_activity_modeling:
#     input:
#         protein_activity_scgenexpr = os.path.join(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_K562_essential-pseudobulk_across_batches-scgenexpr.tsv.gz")
#     output:
#         directory(os.path.join(RESULTS_DIR,"figures","activity_modeling"))
#     shell:
#         """
#         Rscript scripts/figures_activity_modeling.R \
#                     --evaluation_ex_file={input.evaluation_ex} \
#                     --evaluation_genexpr_file={input.evaluation_genexpr} \
#                     --protein_activity_ex_file={input.protein_activity_ex} \
#                     --protein_activity_genexpr_file={input.protein_activity_genexpr} \
#                     --figs_dir={output}
#         """
       