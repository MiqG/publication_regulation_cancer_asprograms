import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","network_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

VIPER_SPLICING_DIR = os.path.join(ROOT,"../../repositories/viper_splicing")


PERT_SPLICING_FILES = {
    "ENCOREKD_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKD_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKO_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'delta_psi-{omic_type}.tsv.gz'),
    "ENASFS": os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','delta_psi-{omic_type}.tsv.gz')
}

PERT_GENEXPR_FILES = {
    "ENCOREKD_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKD_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKO_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'log2_fold_change_tpm.tsv.gz'),
    "ENASFS": os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','log2_fold_change_tpm.tsv.gz'),
    "ReplogleWeissman2022_K562_essential-pseudobulk_across_batches": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_essential-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
    "ReplogleWeissman2022_rpe1-pseudobulk_across_batches": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz"),
    "ReplogleWeissman2022_K562_gwps-pseudobulk_across_batches": os.path.join(PREP_DIR,"pert_transcriptomes","ReplogleWeissman2022_K562_gwps-pseudobulk_across_batches-log2_fold_change_cpm.tsv.gz")
}

PERT_FILES = {
    "EX": PERT_SPLICING_FILES,
    "genexpr": PERT_GENEXPR_FILES
}
OMIC_PERT_DICT = {
    "EX": "EX",
    "genexpr": "genexpr",
    "scgenexpr": "genexpr"
}

METADATA_FILES = [
    os.path.join(PREP_DIR,"metadata","ENCOREKO.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENCOREKD.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz"),
    os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_K562_essential-pseudobulk_across_batches-conditions.tsv.gz"),
    os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_rpe1-pseudobulk_across_batches-conditions.tsv.gz"),
    os.path.join(PREP_DIR,"singlecell","ReplogleWeissman2022_K562_gwps-pseudobulk_across_batches-conditions.tsv.gz")
]

REGULON_DIRS = {
    "genexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-genexpr"),
    "EX": os.path.join(VIPER_SPLICING_DIR,"data","empirical_sf_networks-EX"),
    "scgenexpr": os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-scgenexpr")
}


SHADOWS = ["no"] # bug in viper does not allow shadow correction
N_TAILS = ["two"]

EVENT_TYPES = ["EX"]
OMIC_TYPES = ["genexpr","scgenexpr"] + EVENT_TYPES

EVAL_DATASETS = list(PERT_GENEXPR_FILES.keys())

ALL_DATASETS = [d for d in PERT_GENEXPR_FILES.keys() for o in ["genexpr","scgenexpr"]] + [d for d in PERT_SPLICING_FILES.keys() for o in ["EX"]]
ALL_OMICS = [o for d in PERT_GENEXPR_FILES.keys() for o in ["genexpr","scgenexpr"]] + [o for d in PERT_SPLICING_FILES.keys() for o in ["EX"]]
ALL_SHADOWS = [s for s in SHADOWS for d in ALL_DATASETS]
ALL_N_TAILS = [t for t in N_TAILS for d in ALL_DATASETS]

##### RULES #####
rule all:
    input:
        # prepare evaluation labels
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"),
        
        # evaluate regulons
        ## run
        expand(os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{dataset}-{omic_type}-shadow_{shadow}-{n_tails}_tailed.tsv.gz"), zip, dataset=ALL_DATASETS, omic_type=ALL_OMICS, shadow=ALL_SHADOWS, n_tails=ALL_N_TAILS),
        ## merge
        expand(os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # estimate splicing factor activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_type}.tsv.gz"), zip, dataset=ALL_DATASETS, omic_type=ALL_OMICS),
        
        # train models
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_type}_to_EX","{model_type}","losses.tsv.gz"), omic_type=["genexpr","scgenexpr"], model_type=["fclayer","ewlayer"]),
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_type}_to_EX","{model_type}","weights.pth"), omic_type=["genexpr","scgenexpr"], model_type=["fclayer","ewlayer"]),
        expand(os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_type}_to_EX","{model_type}","common_regulators.tsv.gz"), omic_type=["genexpr","scgenexpr"], model_type=["fclayer","ewlayer"]),
        
        # make predictions
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_{model_type}_and_{omic_type}.tsv.gz"), dataset=PERT_SPLICING_FILES.keys(), omic_type=["genexpr","scgenexpr"], model_type=["fclayer","ewlayer"])
        
# #         # make figures
# #         os.path.join(RESULTS_DIR,"figures","eval_genexpr_vs_splicing")
        
        
rule make_evaluation_labels:
    input:
        metadatas = METADATA_FILES
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"))
    run:
        import pandas as pd
        
        for f in input.metadatas:
            # load
            metadata = pd.read_table(f)
            
            os.makedirs(output.output_dir, exist_ok=True)
            if "ENCORE" in f:
                for cell_line in metadata["cell_line"].unique():
                    # make labels
                    labels = metadata.loc[
                        metadata["cell_line"]==cell_line, ["PERT_GENE","PERT_ENSEMBL"]
                    ].drop_duplicates().copy()
                    labels["PERT_ID"] = metadata["PERT_ENSEMBL"]
                    labels["PERT_TYPE"] = "KNOCKDOWN" if "KD" in os.path.basename(f) else "KNOCKOUT"
                    
                    dataset = "ENCOREKD" if "KD" in os.path.basename(f) else "ENCOREKO"
                    
                    # save
                    labels.dropna().to_csv(os.path.join(output.output_dir,"%s_%s.tsv.gz") % (dataset, cell_line), **SAVE_PARAMS)
            elif "singlecell" in f:
                metadata["PERT_ID"] = metadata["PERT_ENSEMBL"]
                metadata["PERT_TYPE"] = "KNOCKDOWN"
                # "PERT_GENE" columns would be nice to have
                labels = metadata[["PERT_ID","PERT_ENSEMBL","PERT_TYPE"]].drop_duplicates()
                
                # save
                dataset_file = os.path.basename(f).replace("-conditions","")
                labels.dropna().to_csv(os.path.join(output.output_dir,dataset_file), **SAVE_PARAMS)
                
            elif "ENASFS" in f:
                # prepare labels
                metadata["PERT_ID"] = metadata[
                    ["study_accession","cell_line_name","PERT_ENSEMBL"]
                ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
                labels = metadata[["PERT_ID","PERT_GENE","PERT_ENSEMBL","PERT_TYPE"]].drop_duplicates()

                # save
                labels.dropna().to_csv(os.path.join(output.output_dir,"ENASFS.tsv.gz"), **SAVE_PARAMS)
                
                
        print("Done!")
        

rule evaluate_regulons:
    input:
        signature = lambda wildcards: PERT_FILES[OMIC_PERT_DICT[wildcards.omic_type]][wildcards.dataset],
        regulons = lambda wildcards: REGULON_DIRS[wildcards.omic_type],
        eval_labels = os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels","{dataset}.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{dataset}-{omic_type}-shadow_{shadow}-{n_tails}_tailed.tsv.gz")
    params:
        script_dir = SRC_DIR,
        shadow = "{shadow}",
        n_tails = "{n_tails}"
    shell:
        """
        nice Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons} \
                    --eval_labels_file={input.eval_labels} \
                    --output_file={output} \
                    --shadow_correction={params.shadow} \
                    --n_tails={params.n_tails}
        """
        
        
rule combine_evaluations:
    input:
        evaluations = lambda wildcards: [os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{dataset}-{omic_type}-shadow_{shadow}-{n_tails}_tailed.tsv.gz").format(dataset=d, omic_type="{omic_type}", shadow=s, n_tails=n) for s in SHADOWS for n in N_TAILS for d in PERT_FILES[OMIC_PERT_DICT[wildcards.omic_type]].keys()]
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-{omic_type}.tsv.gz")
    params:
        omic_type = "{omic_type}"
    run:
        import pandas as pd
    
        evaluation = pd.concat([pd.read_table(f) for f in input.evaluations])
        evaluation["omic_type"] = params.omic_type
        
        evaluation.to_csv(output[0], **SAVE_PARAMS)
        
        print("Done!")
        
        
rule compute_protein_activity:
    input:
        signature = lambda wildcards:PERT_FILES[OMIC_PERT_DICT[wildcards.omic_type]][wildcards.dataset],
        regulons_path = lambda wildcards: REGULON_DIRS[wildcards.omic_type]
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_type}.tsv.gz")
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
        train_features = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKD_HepG2-{omic_type}.tsv.gz"),
        train_labels = os.path.join(RESULTS_DIR,"files","protein_activity","ENCOREKD_HepG2-EX.tsv.gz"),
        val_features = os.path.join(RESULTS_DIR,"files","protein_activity","ENASFS-{omic_type}.tsv.gz"),
        val_labels = os.path.join(RESULTS_DIR,"files","protein_activity","ENASFS-EX.tsv.gz"),
    output:
        losses = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_type}_to_EX","{model_type}","losses.tsv.gz"),
        weights = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_type}_to_EX","{model_type}","weights.pth"),
        common_regulators = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_type}_to_EX","{model_type}","common_regulators.tsv.gz")
    params:
        model_type = "{model_type}",
        batch_size = 564,
        max_epochs = 50,
        learning_rate = 0.01
    run:
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset
        import lightning as L
        
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
        
        # subset
        common_samples_train = set(train_features.columns).intersection(train_labels.columns)
        common_samples_val = set(val_features.columns).intersection(val_labels.columns)
        common_regulators = set(
            set(train_features.index).intersection(train_labels.index)
        ).intersection(
            set(val_features.index).intersection(val_labels.index)
        )
        train_features = train_features.loc[common_regulators, common_samples_train].copy()
        train_labels = train_labels.loc[common_regulators, common_samples_train].copy()
        val_features = val_features.loc[common_regulators, common_samples_val].copy()
        val_labels = val_labels.loc[common_regulators, common_samples_val].copy()
        
        # prep dataloaders
        X_train = torch.tensor(train_features.fillna(0).T.values, dtype=torch.float32)
        Y_train = torch.tensor(train_labels.fillna(0).T.values, dtype=torch.float32)
        X_val = torch.tensor(val_features.fillna(0).T.values, dtype=torch.float32)
        Y_val = torch.tensor(val_labels.fillna(0).T.values, dtype=torch.float32)
        
        ds_train = TensorDataset(X_train, Y_train)
        ds_val = TensorDataset(X_val, Y_val)
        
        train_loader = DataLoader(ds_train, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(ds_val, batch_size=batch_size, shuffle=False)
        
        # define model
        ## model
        class ElementwiseLinear(nn.Module):
            def __init__(self, input_size: int) -> None:
                super(ElementwiseLinear, self).__init__()

                # w is the learnable weight of this layer module
                self.w = nn.Parameter(torch.rand(input_size), requires_grad=True)

            def forward(self, x: torch.tensor) -> torch.tensor:
                # simple elementwise multiplication
                return self.w * x
        
        class EWlayer(nn.Module):
            def __init__(self, input_size):
                super(EWlayer, self).__init__()
                self.ew = ElementwiseLinear(input_size)
                
            def forward(self, x):
                x = self.ew(x)
                return x

        class FClayer(nn.Module):
            def __init__(self, input_size, output_size):
                super(FClayer, self).__init__()
                self.fc = nn.Linear(input_size, output_size)
                
            def forward(self, x):
                x = self.fc(x)
                return x
            
        ## wraper for training with lightning
        class LitWrapper(L.LightningModule):
            def __init__(self, model, criterion, learning_rate):
                super(LitWrapper, self).__init__()
                self.model = model
                self.criterion = criterion
                self.learning_rate = learning_rate
                self.losses = []
                
            def forward(self, x):
                return self.model(x)

            def training_step(self, batch, batch_idx):
                x, y = batch
                y_pred = self(x)
                loss = self.criterion(y_pred, y)
                self.log('train_loss', loss)
                self.losses.append(
                    {"epoch": self.current_epoch, "batch": batch_idx, "train_loss": float(loss)}
                )
                return loss

            def validation_step(self, batch, batch_idx):
                x, y = batch
                y_pred = self(x)
                loss = self.criterion(y_pred, y)
                self.log('val_loss', loss)
                self.losses.append(
                    {"epoch": self.current_epoch, "batch": batch_idx, "val_loss": float(loss)}
                )
                return loss
            
            def configure_optimizers(self):
                return torch.optim.Adam(self.parameters(), lr=self.learning_rate)
            
        # train model
        criterion = nn.MSELoss()
        
        if model_type=="fclayer":
            model = FClayer(input_size=len(common_regulators), output_size=len(common_regulators))
        
        elif model_type=="ewlayer":
            model = EWlayer(input_size=len(common_regulators))
        
        litwrapper = LitWrapper(model, criterion, learning_rate)
        trainer = L.Trainer(
            max_epochs = max_epochs,
            accelerator = "cpu",
            enable_checkpointing=False, # do not save checkpoints
            logger=False # do not save logs
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
        pd.DataFrame(common_regulators).to_csv(output.common_regulators, header=None, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule predict_sf_activity_from_model:
    input:
        activity = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_type}.tsv.gz"),
        weights = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_type}_to_EX","{model_type}","weights.pth"),
        common_regulators = os.path.join(RESULTS_DIR,"files","model_sf_activity","from_{omic_type}_to_EX","{model_type}","common_regulators.tsv.gz")
    output:
        activity_pred = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-EX_from_model_{model_type}_and_{omic_type}.tsv.gz")
    params:
        model_type = "{model_type}"        
    run:
        import torch
        import torch.nn as nn
        import pandas as pd
        
        # load
        activity = pd.read_table(input.activity, index_col=0)
        weights = torch.load(input.weights)
        common_regulators = list(pd.read_table(input.common_regulators, header=None)[0])
        model_type = params.model_type
        
        # prep data
        activity = pd.merge(
            pd.DataFrame(index=common_regulators),
            activity,
            how="left", left_index=True, right_index=True
        ).fillna(0)
        X = torch.tensor(activity.T.values, dtype=torch.float32)
        
        # load model
        class ElementwiseLinear(nn.Module):
            def __init__(self, input_size: int) -> None:
                super(ElementwiseLinear, self).__init__()

                # w is the learnable weight of this layer module
                self.w = nn.Parameter(torch.rand(input_size), requires_grad=True)

            def forward(self, x: torch.tensor) -> torch.tensor:
                # simple elementwise multiplication
                return self.w * x
        
        class EWlayer(nn.Module):
            def __init__(self, input_size):
                super(EWlayer, self).__init__()
                self.ew = ElementwiseLinear(input_size)
                
            def forward(self, x):
                x = self.ew(x)
                return x

        class FClayer(nn.Module):
            def __init__(self, input_size, output_size):
                super(FClayer, self).__init__()
                self.fc = nn.Linear(input_size, output_size)
                
            def forward(self, x):
                x = self.fc(x)
                return x
            
            
        if model_type=="fclayer":
            model = FClayer(input_size=len(common_regulators), output_size=len(common_regulators))
        
        elif model_type=="ewlayer":
            model = EWlayer(input_size=len(common_regulators))
        
        # make predictions
        model.eval()
        with torch.no_grad():
            Y_hat = model(X)
        
        # prep outputs
        activity_pred = pd.DataFrame(Y_hat.detach().numpy().T, index=common_regulators, columns=activity.columns)
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
       