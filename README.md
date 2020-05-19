# metPropagate

This repository contains code for the metPropagate method, published in npj Genomic Medicine. Citation: 

Emma Graham Linck, Phillip A Richmond, Maja Tarailo-Graovac, Udo Engelke, Leo A.J. Kluijtmans, Karlien L.M Coene, Ron A. Wevers, Wyeth Wasserman, Clara D.M. van Karnebeek Sara Mostafavi. metPropagate: network-guided propagation of metabolomic information for prioritization of metabolic disease genes. npj Genomic Medicine (May 2020).

##Explanation of repository contents

This code in this repo is divided into three sections: 

1. METABOLOMIC_PROCESSING_PIPELINE

    input: 
    - candidate_genes/sampleA_variants.csv: this is a list of genes in which to assess metabolomic enrichment
    - metabolomics/negative_input_intensity_matrix.csv
    - metabolomics/positive_input_intensity_matrix.csv
    - metabolomics/negativepeakTable.csv
    - metabolomics/positivepeakTable.csv
    
    The candidate genes file can be replaced with a list of patient-specific candidate genes. 
    The sample metabolomics files provided here are the output of an XCMS analysis comparing a single patient to a group of controls. The XCMS code used to generate these files is commented out in the first section of the preprocessing.R file. 

    output:
    - neg_mode or pos_mode/linear_raw_intensities.csv (linear scaled input intensity matrix from ESI-/+ mode data only)
    - neg_mode or pos_mode/log_raw_intensities.csv (log scaled input intensity matrix from ESI-/+ mode data only)
    - neg_mode or pos_mode/unnormalized_raw_intensities.csv (input intensity matrix from ESI-/+ mode data only)
    - neg_mode or pos_mode/negative(or positive)_linear_mode_msea_genebig_allpeaks.csv (enrichment of DAMs for all genes in HMDB, given input of linear scaled intensity matrix from ESI-/+ mode data only)
    - neg_mode or pos_mode/number_of_metabolites_matched_to_features.csv (number of possible metabolite IDs of each feature after exact mass matching of input m/z ratios)
    - neg_mode or pos_mode/significant_features_linear.csv (feature IDs of significant features in all samples)
    - output/linear_msea_genebig_allpeaks_combined.csv (enrichment of DAMs for all genes in HMDB, when combining both positive and negative mode DAMs)

2. label_propagation

    This section contains code for the label propagation component of metPropagate. It was adapted from Yuto Yamaguchi (https://github.com/yamaguchiyuto/label_propagation). 
    
    This folder contains the bash script that runs label propagation given inputs from METABOLOMIC_PROCESSING_PIPELINE: run_label_prop_cluster.sh
    
3. integration

    The contents of enrichment_files and wes_files are populated directly from METABOLOMIC_PROCESSING_PIPELINE. 
    
    graph_files contains input STRING network files (manually placed)
        - STRING_graph_file_v11_gene_list_functional_entire_db: STRING database downloaded from stringdb.org
        - STRING_graph_file_v11_gene_to_nodeid_mapping_functional_entire_db: node_id to gene name mapping file. 
    enrichment_files contains input metabolomic enrichment files that are output from METABOLOMIC_PROCESSING_PIPELINE (manually placed)
        - sampleA.csv is a sample enrichment file
    label_files contains per-gene label files based on metabolomic enrichment data used as input to label propagation (automatically generated)
        - STRING_sampleA is an intermediate file. 
        - intermediate_label_file_sampleA is an example of a label file that is input to label propagation 
    wes_files contains candidate gene information
        - sampleA_candidate_genes.csv is an example input file from Exomiser
    LPA_output contains the post-propagation node scores
        - sampleA is example file
        
        
##Quick start
 
 System requirements:
 
 This will get you started with the dummy input files already included in this repo: 
  1. Run *METABOLOMIC_PROCESSING_PIPELINE/run_met_processing.sh*. Output is metabolomic enrichment files in output folder. Note that this file calls .pbs script, and should be modified if you are using a different batch processing system (eg. SLURM).
  2. Move output file named *linear_msea_genebig_allpeaks_combined.csv* in *METABOLOMIC_PROCESSING_PIPELINE/output/* to *integration/enrichment_files*. Rename it *samplename_enrichment_file.csv*. 
  3. Edit the sample_name variable in *label_propagation/run_label_prop_cluster.sh* to be the sample(s) you want to analyze. 
  4. Run *run_label_prop_cluster.sh*. Note that this file calls a .pbs file. Results of label propagation will populate in *integration/results/*. File contains final propagated scores for each candidate gene as "Score.ID" and the original scores for each gene as "Original.Score.ID". 
    
