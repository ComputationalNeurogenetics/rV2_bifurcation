# rV2_bifurcation

This is code and metadata repository for manuscript "Regulation of Tal1 dependent lineage bifurcation in rhombomere 1".

## Processing scRNA data

1. code/e_all_refilter_rescale_complete.Rmd 
(TODO: clean stashed changes, all non-relevant comments and E14/E15 data handling?)
- Datafile is stored in scRNA_data/e_all_rescaled_regress_nuisance.Rds

## Processing scATAC data

1. code/JointFeatureSpace.Rmd (TODO: need to be cleaned and rerun based on notebook in analysis/ to generate clean Rmd and html)
- Joined featurespace is stored in analysis/JointFeatureSpace.271021.Rds (equivalent file is joint.clade.peaks.h8.281021.Rds, TODO: clean extra)

2. code/E12R1_clade_peaks_joint_space_271021.Rmd

3. code/e12r_counts_to_downstream.240322.Rmd
- Datafile for downstream processing is stored in scATAC_data/E12_R1_DownstreamReady_nmm_.240322.Rds

## rV2-lineage similarity pattern between E12.5 and E13.5

1. code/gaba_cross_sample_cluster_comparison_monocle3.Rmd and code/glut_cross_sample_cluster_comparison_monocle3.Rmd

## VIA pseudotime analysis

## TF-footprinting

## Cut & Run data processing
