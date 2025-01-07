# rV2_bifurcation

This is code and metadata repository for manuscript "Regulation of Tal1 dependent lineage bifurcation in rhombomere 1".

## Processing scRNAseq
1. code/e_all_refilter_rescale_complete.Rmd 
    - Results in (scRNA_data/e_all_rescaled_regress_nuisance.Rds, but not in Github)

## Processing scATACseq
1. code/E12R1_clade_peaks_joint_space_271021.Rmd
    - Featurespace is stored in analysis/JointFeatureSpace.271021.Rds
2. code/e12r_counts_to_downstream.240322.Rmd
    - Results in (scATAC_data/E12_R1_DownstreamReady_nmm_.240322.Rds, but not in Github)

## rV2-lineage subsetting
1. label transfer rna_integration_atac_label_anno_transfer.Rmd
2. rV2-lineage labeling rV2_lineage_labeling.Rmd
    FindSubCluster for aggregated GA/GL half-lineages, initial subset based on FindSubCluster-generated community labels
    results in scATAC_data/nmm_rV2_subset_relabeled_110522.Rds, but not in Github
3.nmm_rV2_subset_relabeled_110522.Rds => nmm_rV2_subset_relabeled_031023_links.qs with updated links information

## VIA pseudotime analysis for scRNAseq cells
1. code/VIA_pseudotime_scRNA.Rmd
2. code/E12_rV2_scRNA_200923.ipynb

## VIA pseudotime analysis for scATACseq cells
1. code/VIA_preparation.Rmd
2. code/rV2_VIA_pseudotime.ipynb

## TF-footprinting


## RNAscope image segmentation and analysis
### StarDist based segmentation
1.Creation of training images for StarDist model building
- 10 images covering majority of cell shapes and sizes in the dissections of E12.5 R1 were manually painted based on DAPI channel information to identify each cell. 
2. Training images were used to train Startdist model as per /code/RNAscope/Stardist_training.ipynb
3. DAPI channel of each actual analysis images was then used to segment images and discover cell boundaries as per /code/RNAscope/Prediction.ipynb
4. ROIs of each image were exported as ImageJ ROIS.zip files

### Downstream analysis of intensities of fluorescence channels in ROIs
1. ROIs of each image were read into R with intensities per channel per image as per /code/RNAscope/SegmentationAnalysis_v2.Rmd
2. Futher downstream analysis to calculate intensity changes of each channel as a function of cell migration alogn radial axis (away from ventricular zone) were calculated as per /code/RNAscope/SegmentationAnalysis_v2.Rmd

## Cut & Run data processing
1. nf-core/cutandrun pipeline with run report (code/CT/pipeline_report.txt), including launch parameters



