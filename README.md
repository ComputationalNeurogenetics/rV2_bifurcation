# rV2_bifurcation

This is code and metadata repository for manuscript "Regulation of Tal1 dependent lineage bifurcation in rhombomere 1".

## Processing scRNA data

1. code/e_all_refilter_rescale_complete.Rmd 
(TODO: clean stashed changes, all non-relevant comments and E14/E15 data handling?)
- Datafile is stored in scRNA_data/e_all_rescaled_regress_nuisance.Rds

## rV2-lineage similarity pattern between E12.5 and E13.5

1. code/gaba_cross_sample_cluster_comparison_monocle3.Rmd and code/glut_cross_sample_cluster_comparison_monocle3.Rmd

## Processing scATAC data
Mistä RNA data tuli
1. code/JointFeatureSpace.Rmd (TODO: need to be cleaned and rerun based on notebook in analysis/ to generate clean Rmd and html)
- Joined featurespace is stored in analysis/JointFeatureSpace.271021.Rds (equivalent file is joint.clade.peaks.h8.281021.Rds, TODO: clean extra)

2. code/E12R1_clade_peaks_joint_space_271021.Rmd

3. code/e12r_counts_to_downstream.240322.Rmd
- Datafile for downstream processing is stored in scATAC_data/E12_R1_DownstreamReady_nmm_.240322.Rds

## rV2-lineage subsetting

0. label transfer https://github.com/laasov/Gradu/blob/master/code/rna_integration_atac_label_anno_transfer.Rmd (transferred labels not used) labels from document cited in code

1. rV2-lineage labeling https://github.com/laasov/Gradu/blob/master/code/rV2_lineage_labeling.Rmd
    FindSubCluster for aggregated GA/GL half-lineages, initial subset based on FindSubCluster-generated community labels
    data object nmm_rV2_subset_relabeled_110522.Rds

2.


## VIA pseudotime analysis

1. code/VIA_preparation.Rmd
2. code/rV2_VIA_pseudotime.ipynb
3. VIA_pseudotime_features.Rmd

## TF-footprinting

## Cut & Run data processing

1. nf-core/cutandrun quality filtering and alignment steps
2. MACS2 (https://github.com/macs3-project/MACS) with parameters
macs2 callpeak [target] [control] -f BAMPE --nomodel -g mm --gsize 1.87e9 --nolambda -q 0.1 --max-gap 10 --tsize 8

 IgG treated samples were used as controls. Number of replicates for consensus peaks: n(Gata2)=4, n(Gata3)=3, n(H3K4me3)=6


