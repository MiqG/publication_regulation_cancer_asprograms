import os

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
        "genexpr": os.path.join(RAW_DIR,"viper_splicing_intermediate_files","genexpr_tpm","CCLE.tsv.gz"), 
        "EX": os.path.join(RAW_DIR,"viper_splicing_intermediate_files","event_psi","CCLE-EX.tsv.gz")
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

OMIC_TYPES = ["bulkgenexpr","scgenexpr","bulkscgenexpr","EX"]
OMIC_GENEXPR_REGULONS = ["bulkgenexpr","scgenexpr","bulkscgenexpr"]
MODEL_ARCHS = ["fclayer","ewlayer"]
K_CROSS_VALIDATION = 5

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
        [os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights-{k}.pth").format(k=k, omic_regulon=o, model_type=m) for o in OMIC_GENEXPR_REGULONS for m in MODEL_ARCHS for k in range(K_CROSS_VALIDATION)],
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","input_regulators.tsv.gz"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_ARCHS),
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","output_regulators.tsv.gz"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_ARCHS),
        os.path.join(RESULTS_DIR,"files","model_sf_activity","losses-merged.tsv.gz"),
        
        # make figures
        os.path.join(RESULTS_DIR,"figures","activity_modeling")
        
        
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
        weights = [os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights-{k}.pth").format(k=k, omic_regulon="{omic_regulon}", model_type="{model_type}") for k in range(K_CROSS_VALIDATION)],
        input_regulators = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","input_regulators.tsv.gz"),
        output_regulators = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","output_regulators.tsv.gz")
    params:
        model_type = "{model_type}",
        batch_size = 128,
        max_epochs = 20,
        learning_rate = 0.01,
        train_split = 0.85,
        k_cross_validation = K_CROSS_VALIDATION
    threads: 4
    run:
        import numpy as np
        import pandas as pd
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset
        import lightning as L
        from vipersp.model import EWlayer, FClayer, LitWrapper
        from sklearn.model_selection import KFold
        
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
        k_cross_validation = params.k_cross_validation
        
        # subset
        ## samples
        avail_samples = list(set(features.columns).intersection(labels.columns))
        ## regulators
        if model_type=="fclayer":
            input_regulators = list(features.index)
            output_regulators = list(labels.index)
        elif model_type=="ewlayer":
            common_regulators = list(set(features.index).intersection(labels.index))
            input_regulators = common_regulators
            output_regulators = common_regulators
        
        kf = KFold(n_splits=k_cross_validation, shuffle=True, random_state=1234)
        
        all_losses = []
        for k, (train_index, test_index) in enumerate(kf.split(features.T)):
            # split
            train_samples = np.array(avail_samples)[train_index]
            val_samples = np.array(avail_samples)[test_index]
            
            train_features = features.loc[input_regulators, train_samples].copy()
            train_labels = labels.loc[output_regulators, train_samples].copy()
            val_features = features.loc[input_regulators, val_samples].copy()
            val_labels = labels.loc[output_regulators, val_samples].copy()

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
                model = FClayer(input_size=len(input_regulators), output_size=len(output_regulators))

            elif model_type=="ewlayer":
                model = EWlayer(input_size=len(input_regulators))        

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
            losses["k_cross_validation"] = k

            # save
            ## loss
            all_losses.append(losses)
            ## model weights        
            torch.save(model.state_dict(), output.weights[k])
        
        # save
        ## loss
        all_losses = pd.concat(all_losses)
        all_losses.to_csv(output.losses, **SAVE_PARAMS)
        ## common regulators
        pd.DataFrame(input_regulators).to_csv(output.input_regulators, header=None, **SAVE_PARAMS)
        pd.DataFrame(output_regulators).to_csv(output.output_regulators, header=None, **SAVE_PARAMS)
        
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
        
        
rule figures_activity_modeling:
    input:
        losses = os.path.join(RESULTS_DIR,"files","model_sf_activity","losses-merged.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","activity_modeling"))
    shell:
        """
        Rscript scripts/figures_activity_modeling.R \
                    --losses_file={input.losses} \
                    --figs_dir={output}
        """