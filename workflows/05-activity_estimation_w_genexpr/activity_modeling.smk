import os
import copy

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","network_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

VIPER_SPLICING_DIR = os.path.join(ROOT,"../../repositories/viper_splicing")

DATASETS = {
    "CCLE": {
        "genexpr": os.path.join(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz"), 
        "EX": os.path.join(PREP_DIR,"event_psi","CCLE-EX.tsv.gz")
    }
}

PERT_SPLICING_FILES = {
    "ENCOREKD_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'delta_psi-EX.tsv.gz'),
    "ENCOREKD_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'delta_psi-EX.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'delta_psi-EX.tsv.gz'),
    "ENCOREKO_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'delta_psi-EX.tsv.gz'),
    "ENASFS": os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','delta_psi-EX.tsv.gz'),
    "Rogalska2024": os.path.join(PREP_DIR,'delta_psi','Rogalska2024-EX.tsv.gz')
}

PERT_GENEXPR_FILES = {
    "ENCOREKD_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKD_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKO_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'log2_fold_change_tpm.tsv.gz'),
    "ENASFS": os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','log2_fold_change_tpm.tsv.gz'),
    "Rogalska2024": os.path.join(PREP_DIR,'log2_fold_change_tpm','Rogalska2024-genexpr_tpm.tsv.gz'),
    "ReplogleWeissman2022_K562_essential-pseudobulk_across_batches": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_essential-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
    "ReplogleWeissman2022_rpe1-pseudobulk_across_batches": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
    "ReplogleWeissman2022_K562_gwps-pseudobulk_across_batches": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_gwps-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz")
}

PERT_EVAL_FILES = {
    "EX": PERT_SPLICING_FILES,
    "genexpr": PERT_GENEXPR_FILES
}

# add CCLE
PERT_FILES = copy.deepcopy(PERT_EVAL_FILES)
PERT_FILES["EX"]["CCLE"] = os.path.join(RESULTS_DIR,"files","signatures","CCLE-EX.tsv.gz")
PERT_FILES["genexpr"]["CCLE"] = os.path.join(RESULTS_DIR,"files","signatures","CCLE-genexpr.tsv.gz")

OMIC_PERT_DICT = {
    "EX": "EX",
    "genexpr": "genexpr",
    "scgenexpr": "genexpr"
}

METADATA_FILES = [
    os.path.join(PREP_DIR,"metadata","ENCOREKO.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENCOREKD.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","Rogalska2024.tsv.gz"),
    os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_K562_essential-pseudobulk_across_batches-conditions.tsv.gz"),
    os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-conditions.tsv.gz"),
    os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_K562_gwps-pseudobulk_across_batches-conditions.tsv.gz")
]

REGULON_SETS = [
    "experimentally_derived_regulons_pruned-genexpr",
    "experimentally_derived_regulons_pruned-scgenexpr",
    "experimentally_derived_regulons_pruned_w_viper_networks-EX"
]

METHODS_ACTIVITY = ["viper","correlation_pearson","correlation_spearman","gsea"]

EVENT_TYPES = ["EX"]
OMIC_TYPES = ["genexpr","scgenexpr"] + EVENT_TYPES
MODEL_TYPES = ["fclayer","ewlayer"]
OMIC_GENEXPR_REGULONS = ["genexpr","scgenexpr"]

EVAL_DATASETS = list(PERT_GENEXPR_FILES.keys())

# ALL_DATASETS = [d for d in PERT_GENEXPR_FILES.keys() for o in OMIC_GENEXPR_REGULONS] + [d for d in PERT_SPLICING_FILES.keys() for o in ["EX"]]
# ALL_OMICS = [o for d in PERT_GENEXPR_FILES.keys() for o in OMIC_GENEXPR_REGULONS] + [o for d in PERT_SPLICING_FILES.keys() for o in ["EX"]]

##### RULES #####
rule all:
    input:
        # signature CCLE
        expand(os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_signature}.tsv.gz"), dataset=["CCLE"], omic_signature=["EX","genexpr"]),

        # estimate splicing factor activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz"), zip, dataset=ALL_DATASETS, omic_regulon=ALL_OMICS),
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz"), dataset=["CCLE"], omic_regulon=OMIC_TYPES),
        
        # train models
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","losses.tsv.gz"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights.pth"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","genes_train.tsv.gz"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","regulators_train.tsv.gz"), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_TYPES),
        
        # estimate SF activities with trained models
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_{model_type}_and_{omic_regulon}.tsv.gz"), dataset=PERT_SPLICING_FILES.keys(), omic_regulon=OMIC_GENEXPR_REGULONS, model_type=MODEL_TYPES),

        # make figures
        # os.path.join(RESULTS_DIR,"figures","activity_modeling")
        
        
rule compute_signature_within:
    input:
        data = lambda wildcards: DATASETS[wildcards.dataset][wildcards.omic_signature]
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_signature}.tsv.gz")
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
        signature = lambda wildcards: PERT_FILES[OMIC_PERT_DICT[wildcards.omic_regulon]][wildcards.dataset],
        regulons_path = lambda wildcards: REGULON_DIRS[wildcards.omic_regulon]
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
        """
        
        
rule train_model_to_predict_sf_activity:
    input:
        train_features = os.path.join(RESULTS_DIR,"files","protein_activity","CCLE-{omic_regulon}.tsv.gz"),
        train_labels = os.path.join(RESULTS_DIR,"files","protein_activity","CCLE-EX.tsv.gz"),
        # train_features = os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'log2_fold_change_tpm.tsv.gz'),
        # train_labels = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKO_HepG2-EX.tsv.gz"),
        # train_features = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKO_HepG2-{omic_regulon}.tsv.gz"),
        # train_labels = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKO_HepG2-EX.tsv.gz"),
        # val_features = os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','log2_fold_change_tpm.tsv.gz'),
        val_features = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKO_HepG2-{omic_regulon}.tsv.gz"),
        val_labels = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKO_HepG2-EX.tsv.gz"),
    output:
        losses = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","losses.tsv.gz"),
        weights = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights.pth"),
        genes_train = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","genes_train.tsv.gz"),
        regulators_train = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","regulators_train.tsv.gz")
    params:
        model_type = "{model_type}",
        batch_size = 64,
        max_epochs = 20,
        learning_rate = 0.01
    threads: 5
    run:
        import pandas as pd
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset
        import lightning as L
        from vipersp.model import EWlayer, FClayer, LitWrapper
        
        # load data
        train_features = pd.read_table(input.train_features, index_col=0)
        train_labels = pd.read_table(input.train_labels, index_col=0)
        val_features = pd.read_table(input.val_features, index_col=0)
        val_labels = pd.read_table(input.val_labels, index_col=0)
        
        # unpack params
        model_type = params.model_type
        batch_size = params.batch_size
        max_epochs = params.max_epochs
        learning_rate = params.learning_rate
        num_workers = threads
        
        # subset
        ## training
        common_samples_train = set(train_features.columns).intersection(train_labels.columns)
        common_index = set(train_features.index).intersection(train_labels.index)
        genes_train = common_index
        regulators_train = common_index
        train_features = train_features.loc[genes_train, common_samples_train].copy()
        train_labels = train_labels.loc[regulators_train, common_samples_train].copy()
        
        ## validation (put everything in the same order as in training)
        common_samples_val = set(val_features.columns).intersection(val_labels.columns)
        val_features = pd.merge(
            pd.DataFrame(index=genes_train),
            val_features[common_samples_val],
            how="left", left_index=True, right_index=True
        )
        val_labels = pd.merge(
            pd.DataFrame(index=regulators_train),
            val_labels[common_samples_val],
            how="left", left_index=True, right_index=True
        )  
        
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
            model = FClayer(input_size=len(genes_train), output_size=len(regulators_train))
        
        elif model_type=="ewlayer":
            model = EWlayer(input_size=len(genes_train))        
            
        
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
        losses = losses.melt(id_vars="epoch", value_vars=["val_loss","train_loss"], var_name="loss_type", value_name="loss")
        
        # save
        ## loss
        losses.to_csv(output.losses, **SAVE_PARAMS)
        ## model weights        
        torch.save(model.state_dict(), output.weights)
        ## common regulators
        pd.DataFrame(genes_train).to_csv(output.genes_train, header=None, **SAVE_PARAMS)
        pd.DataFrame(regulators_train).to_csv(output.regulators_train, header=None, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule predict_sf_activity_from_model:
    input:
        activity = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_regulon}.tsv.gz"),
        weights = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","weights.pth"),
        genes_train = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","genes_train.tsv.gz"),
        regulators_train = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_regulon}_to_EX","{model_type}","regulators_train.tsv.gz")
    output:
        activity_pred = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_{model_type}_and_{omic_regulon}.tsv.gz")
    params:
        model_type = "{model_type}"        
    run:
        import torch
        import pandas as pd
        from vipersp.model import EWlayer, FClayer
        
        # load
        activity = pd.read_table(input.activity, index_col=0)
        weights = torch.load(input.weights)
        genes_train = list(pd.read_table(input.genes_train, header=None)[0])
        regulators_train = list(pd.read_table(input.regulators_train, header=None)[0])
        model_type = params.model_type
        
        # prep data
        activity = pd.merge(
            pd.DataFrame(index=genes_train),
            activity,
            how="left", left_index=True, right_index=True
        )
        X = torch.tensor(activity.fillna(0).T.values, dtype=torch.float32)
        
        if model_type=="fclayer":
            model = FClayer(input_size=len(genes_train), output_size=len(regulators_train))
        
        elif model_type=="ewlayer":
            model = EWlayer(input_size=len(genes_train))        
            
        model.load_state_dict(weights)        
        
        # make predictions
        model.eval()
        with torch.no_grad():
            Y_hat = model(X)
        
        # prep outputs
        activity_pred = pd.DataFrame(Y_hat.detach().numpy().T, index=regulators_train, columns=activity.columns)
        activity_pred.index.name = "regulator"
        
        # save
        activity_pred.reset_index().to_csv(output.activity_pred, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule figures_network_evaluation:
    input:
        evaluation_ex = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-EX.tsv.gz"),
        evaluation_genexpr = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-genexpr.tsv.gz"),
        evaluation_scgenexpr = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-scgenexpr.tsv.gz"),
        protein_activity_ex = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKD_K562-EX.tsv.gz"),
        protein_activity_genexpr = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKD_K562-genexpr.tsv.gz"),
        protein_activity_scgenexpr = os.path.join(RESULTS_DIR,"files","protein_activity","ReplogleWeissman2022_K562_essential-pseudobulk_across_batches-scgenexpr.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","eval_genexpr_vs_splicing"))
    shell:
        """
        Rscript scripts/figures_eval_genexpr_vs_splicing.R \
                    --evaluation_ex_file={input.evaluation_ex} \
                    --evaluation_genexpr_file={input.evaluation_genexpr} \
                    --protein_activity_ex_file={input.protein_activity_ex} \
                    --protein_activity_genexpr_file={input.protein_activity_genexpr} \
                    --figs_dir={output}
        """
       