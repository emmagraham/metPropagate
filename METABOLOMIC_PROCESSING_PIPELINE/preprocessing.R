#FULL PIPELINE for mzData

#LOAD LIBRARIES
unlink(".RData")
options(BIOCINSTALLER_ONLINE_DCF = FALSE)

#change package locations below
.libPaths("/zfs3/users/grahamemma9/grahamemma9/R/3.5/lib64/R/library")
library(norm)
#library(xcms, lib.loc = "/zfs3/users/grahamemma9/grahamemma9/R/x86_64-pc-linux-gnu-library/3.5")
#library(CAMERA, lib.loc = "/zfs3/users/grahamemma9/grahamemma9/R/3.5/lib64/R/library")
#library(XML, lib.loc = "/zfs3/users/grahamemma9/grahamemma9/R/3.5/lib64/R/library")
library(reshape2, lib.loc = "/zfs3/users/grahamemma9/grahamemma9/R/3.5/lib64/R/library")
library(tidyverse, lib.loc = "/zfs3/users/grahamemma9/grahamemma9/R/3.5/lib64/R/library")
args=(commandArgs(TRUE))
sample_name <-as.character(args[[1]])



#This section obtains a list of all control and patient files, which is stored in the variables: control_files and patient_files.
#All of the code in this script runs twice - once to process all of the metabolomic data in positive ion mode, and once to process
#it in negative ion mode (hence the for loop). The results are combined at the end for a combined metabolic enrichment test.
 for (mode in c("positive", "negative")){
#   if(mode == "positive"){
#     control_files <- list.files(path = "./metabolomics/Controls/Controls_mzData/Files_pos", recursive = TRUE, full.names = TRUE)
#     patient_files <- list.files(path = paste0(c("metabolomics/Patients/", sample_name, "/Files_pos"), collapse = ""), recursive = TRUE, full.names = TRUE)
#   } else{
#     control_files <- list.files(path = "./metabolomics/Controls/Controls_mzData/Files_neg", recursive = TRUE, full.names = TRUE)
#     patient_files <- list.files(path = paste0(c("metabolomics/Patients/", sample_name, "/Files_neg"), collapse = ""), recursive = TRUE, full.names = TRUE)
#   }
#   print(control_files)
#   print(patient_files)
#   
#   our_files <-  c(control_files, patient_files)
#   #XCMS------------------------------------
#   
#   #create data frame of file name, and sample grouping (Control or patient)
#   meta_data <- data.frame(sample_name = sub(basename(our_files), pattern = ".mzData",
#   																					replacement = "", fixed = TRUE),
#   												sample_group = c(rep("Control", 15), sample_name),
#   												stringsAsFactors = FALSE) 
#   #check to see if RDS file containing processed data is available. If yes, load data and proceed to filtering step. 
#   #if not, proceed with processing. 
#   if(file.exists(paste0(c("./metabolomics/Patients/",sample_name, "/",sample_name,"_", mode, "_xsaFA.RDS"), collapse = "")) == TRUE){
#     print("yes")
#     xsaFA <- readRDS(paste0(c("./metabolomics/Patients/",sample_name, "/",sample_name,"_", mode, "_xsaFA.RDS"), collapse = ""))
#     print("got past load RDS")
#   }else{
#     print("no")
#     
#     #XCMS - Peak picking and RT adjustment. Code is from XCMS package tutorial. Parameters are optimized with IPO package.
#     raw_data <- readMSData(files = our_files, pdata = new("NAnnotatedDataFrame", meta_data),
#   											 mode = "onDisk") 
#     cwp <- CentWaveParam(peakwidth = c(3, 80), 
#     										 noise = 1000,
#     										 ppm = 15,
#     										 snthresh = 10,
#     										 mzdiff = 0.0045,
#     										 prefilter = c(3, 100))
#   xdata <- findChromPeaks(raw_data, param = cwp)
#   pdp <- PeakDensityParam(sampleGroups = xdata$sample_group)
#   xdata <- groupChromPeaks(xdata, param = pdp)
#   xdata <- adjustRtime(xdata, param = ObiwarpParam(gapInit = 1.2736,
#   																								 gapExtend = 3.3336))
#   pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
#   												minFraction = 1, 
#   												bw = 0.25,
#   												minSamples = 1)
#   xdata <- groupChromPeaks(xdata, param = pdp)
#   any(is.na(chromPeaks(xdata)))
#   xdata <- fillChromPeaks(xdata)
#   
#   intensity_matrix <- featureValues(xdata, method = "medret", value = "into")
#   feature_def <- featureDefinitions(xdata)
#   write.csv(intensity_matrix, paste0(c("./metabolomics/Patients/", sample_name, "/", mode, "_intensity_matrix_April2019_medret_into.csv"), collapse = ""))
#   write.csv(feature_def, paste0(c("./metabolomics/Patients/", sample_name, "/", mode, "_feature_definitions_April2019_medret_into.csv"), collapse = ""))
# 
#   print("got to xcmsSet")
#   
#   #CAMERA-------------------------------------
#   #This section annoates isotopes and adducts of each feature using CAMERA package. Parameters are default. 
#   
#   #Create an xsAnnotate object
#   xcmsSetObject <- as(xdata, "xcmsSet")
#   xcmsSetObject <- fillPeaks(xcmsSetObject)
#   xsa <- xsAnnotate(xcmsSetObject, sample = c(1:length(our_files)), polarity = mode)
#   #Group after RT value of the xcms grouped peak
#   xsaF <- groupFWHM(xsa, perfwhm=0.6)
#   #Verify grouping
#   xsaC <- groupCorr(xsaF)
#   #Annotate isotopes, could be done before groupCorr
#   xsaFI <- findIsotopes(xsaC)
#   #Annotate adducts
#   xsaFA <- findAdducts(xsaFI, polarity=mode)
#   saveRDS(xsaFA, file = paste0(c("./metabolomics/Patients/",sample_name, "/",sample_name,"_", mode, "_xsaFA.RDS"), collapse = ""))
#   }
#   #Get final peaktable and store
#   print("xsaFA")
#   peaks <- getPeaklist(xsaFA, intval = "intb") %>% tbl_df()
#   
#   colnames(peaks) <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", 
#                                                        "Control", sample_name, unlist(lapply(control_files, function(x){substring(x, nchar(x)-18, nchar(x))})), unlist(lapply(patient_files, function(x){substring(x, nchar(x)-18, nchar(x))})), 
#                                                       "isotopes",
#                                                       "adduct", "pcgroup")
#   peakTable_columns <- peaks %>% 
#     dplyr::select(mz, mzmin, mzmax, rt, rtmin, rtmax, npeaks, Control, sample_name, isotopes, adduct, pcgroup)
#   intensity_matrix <- read.csv(paste0(c("./metabolomics/Patients/", sample_name, "/", mode, "_intensity_matrix_April2019_medret_into.csv"), collapse = ""),
#                                stringsAsFactors = FALSE, 
#                                row.names = 1)
#   peakTable <- cbind(peakTable_columns, intensity_matrix[, c(unlist(lapply(control_files, function(x){substring(x, nchar(x)-18, nchar(x))})), unlist(lapply(patient_files, function(x){substring(x, nchar(x)-18, nchar(x))})))])
#   peakTable <- peakTable[is.finite(rowSums(peakTable[, c(unlist(lapply(control_files, function(x){substring(x, nchar(x)-18, nchar(x))})),
#                                        unlist(lapply(patient_files, function(x){substring(x, nchar(x)-18, nchar(x))})))])),]
#   print("done CAMERA")
#   
#   #STATS--------------------------
#   
#   #set working directory to patient-specific sub-folder in Patients folder
#   setwd(paste0(c("./metabolomics/Patients/", sample_name), collapse = ""))
#   
   if (mode == "positive"){
     if(dir.exists("./output/pos_mode")){
       setwd("./output/pos_mode")
     } else{
       dir.create("./output/pos_mode")
       setwd("./output/pos_mode")
     }
   } else{
     if(dir.exists("./output/neg_mode")){
       setwd("./output/neg_mode")
     } else{
       dir.create("output/neg_mode")
       setwd("./output/neg_mode")
     }
   }
  
  #get data frame of intensities from XCMS DIFFREPORT
  intensity_matrix <- read.csv(paste0(c("../../metabolomics/", mode, "_input_intensity_matrix.csv"), collapse = ""),
                                stringsAsFactors = FALSE, 
                               row.names = 1)
  peakTable <- read.csv(paste0(c("../../metabolomics/", mode, "peakTable.csv"), collapse = ""))
  write.csv(intensity_matrix, "./unnormalized_raw_intensities.csv")
  unnormalized_data <- intensity_matrix[is.finite(rowSums(log2(intensity_matrix))),]
  log_data <- log2(unnormalized_data)[is.finite(rowSums(log2(unnormalized_data))),]
  write.csv(log_data, "./log_raw_intensities.csv")
  
  #linear baseline normalization from Bolstad et al, 2003
  linear.baseline <- apply(log_data,1, function(x){median(x[complete.cases(x)])}) 
  baseline.mean <- mean(linear.baseline, na.rm = TRUE)
  sample.means <- colMeans(log_data, na.rm = TRUE) 
  linear.factor <- baseline.mean/sample.means 
  linear.normalized.data <- t(t(log_data)*linear.factor)
  write.csv(linear.normalized.data, "./linear_raw_intensities.csv")

  print("done with normalization")

  #this function will generate a table with a mz, rt, isotope and adduct column. The input is the 
  #intensity matrix,  an index of peaks to remove, and the table with adduct and isotope annotation.
  addColumns <- function(table_wo_bad_peaks, index, table_with_annotation, remove_peaks){
    if (remove_peaks){
      peaTable_wo_bad_peaks <- table_with_annotation[-c(index),]
    } else{
      peaTable_wo_bad_peaks <- table_with_annotation
    }
    data.frame(table_wo_bad_peaks) %>% dplyr::mutate(mz = peaTable_wo_bad_peaks[["mz"]], 
                                                     rt = peaTable_wo_bad_peaks[["rt"]], 
                                                     isotopes = peaTable_wo_bad_peaks[["isotopes"]],
                                                     adduct = peaTable_wo_bad_peaks[["adduct"]])
  }
  
  #add columns to data with all peaks
 
  linear_normalized_w_columns <- addColumns(linear.normalized.data, badPeaksIndex, peakTable, FALSE)
  print("linear_normalized_w_columns")
  unnormalized_w_columns <- addColumns(unnormalized_data, badPeaksIndex, peakTable, FALSE)
  print("unnormalized_w_columns")
  
  #Isolate annotated base mass from adduct column. If a feature could be multiple base masses, all are included. 
  #If no adduct is annotated to a feature, then the unannotated peak mz is used in its place. 
  addMasses <- function(data_wo_bad_peaks_annotated){
    i <- 1
    mass_vector <- c()
    for (entry in data_wo_bad_peaks_annotated$adduct){
      ids_final <- c()
      if (entry != ""){
        split1 <- unlist(strsplit(entry, split = " "))
        for (x in 1:(length(split1)/2)){
          ids_final <- c(ids_final, split1[x*2])
        }
        if (length(ids_final) > 1){
          mass_vector[i] <- paste0(ids_final, collapse = ";")
          i <- i + 1
        }else{
          mass_vector[i] <- ids_final
          i <- i + 1
        }
      } else{
        mass_vector[i] <- data_wo_bad_peaks_annotated$mz[i]
        i <- i + 1
      }
    }
    return(cbind(data_wo_bad_peaks_annotated, mass = mass_vector))
  }
  
  
  linear_normalized_w_columns_w_masses <- addMasses(linear_normalized_w_columns)
  print("linear_normalized_w_columns_w_masses")
  unnormalized_w_columns_w_masses <- addMasses(unnormalized_w_columns)
  print("unnormalized_w_columns_w_masses")
  rm(linear_normalized_w_columns)
  rm(unnormalized_w_columns)
  #this function removes isotopes that are not base isotopes ([M]+)
  removeIsotopes <- function(data_wo_bad_peaks_annotated_w_masses){
    i <- 1
    discard_isotope_list <- c()
    for (entry in data_wo_bad_peaks_annotated_w_masses$isotopes){
      isotope_l <- grepl("\\[M\\+\\d+\\]", entry)
      if (isotope_l == TRUE){
        discard_isotope_list <- c(discard_isotope_list, i)
        i <- i+ 1
      } else{
        i <- i +1
      }
    }
    return(data_wo_bad_peaks_annotated_w_masses[-c(discard_isotope_list),])
  }
  
  
  linear_normalized_w_columns_w_masses_wo_iso <- removeIsotopes(linear_normalized_w_columns_w_masses)
  print("linear_normalized_w_columns_w_masses_wo_iso")

  unnormalized_w_columns_w_masses_wo_iso <- removeIsotopes(unnormalized_w_columns_w_masses)
  print("unnormalized_w_columns_w_masses_wo_iso")
  
  print("finished data cleaning")
  
  rm(linear_normalized_w_columns_w_masses)
  rm(unnormalized_w_columns_w_masses)
  #ISOLATE FEATURES THAT MAP TO COMPOUNDS
  #read in csv file containing HMDB database
  
  HMDB <- read.csv("../../HMDB_db_ver4.csv", header = TRUE)
  
  #isolate control and patient names
  control_names <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10")
  patient_names <- "IEM1"

  linear_data <- linear_normalized_w_columns_w_masses_wo_iso
  unn_data <- unnormalized_w_columns_w_masses_wo_iso
  rm(linear_normalized_w_columns_w_masses_wo_iso)
  rm(unnormalized_w_columns_w_masses_wo_iso)
  #function that checks whether a given mass has a metabolite with mass within 15ppm in HMDB.
  matchCompounds <- function(imp_peaks_df, db){
    match_or_not <- unlist(lapply(imp_peaks_df$mass, 
                                  function(s){
                                    masses <- as.numeric(unlist(strsplit(as.character(s), split = ";")))
                                    length <- unlist(lapply(masses, function(x){
                                      mass_range <- c(x-(x*(15/1e6)), x+(x*(15/1e6)))
                                      filtered_db <- db %>% filter(monisotopic_molecular_weight >= mass_range[1], monisotopic_molecular_weight <= mass_range[2])
                                      nrow(filtered_db)}))
                                    
                                    sum(length) > 0
                                  }) )
    imp_peaks_df[match_or_not,]
  }
  
  
  linear_compounds_all_peaks1 <- matchCompounds(linear_data, HMDB)
  print("linear_compounds_all_peaks1")
  unn_compounds_all_peaks1 <- matchCompounds(unn_data, HMDB)
  print("unn_compounds_all_peaks1")
  
  print("end preparing metabolite lists")
  rm(linear_data)
  #ISOLATING DIFFERENTIALLY ABUNDANT FEATURES----------------
  
  #this section isolates features that have a z-score greater than 2. Z-score is built using 15 control samples.
  #Change the column that is indexed (16 in code below) depending on how many controls you have. So if you have 20 controls, you would isolate column 21. 
  linear_data1 <- t(linear_compounds_all_peaks1[,colnames(linear_compounds_all_peaks1) %in% c(control_names, patient_names)])
  scaled <- apply(linear_data1, 2, function(x){scale(x, center = TRUE, scale = TRUE)})
  str(scaled)
  special <- abs(scaled[11, ]) >= 2
  print(dim(linear_compounds_all_peaks1))
  print(length(special))
  linear_compounds_all_peaks1$zscore <- scaled[11,]
  linear_sig_matrix_all_peaks <- linear_compounds_all_peaks1[special,]
  write.csv(linear_sig_matrix_all_peaks, "./significant_features_linear.csv")
  
  print("done with stats")
  #HMDB lookup-------------------------------
  
  #For all signficant features, this function identifies the metabolites annotated to each feature.
  #It also identifies the number of metabolites that map to a particular feature, as well as the maximum feature-based 
  #z-score associated
  #with a particular metabolite. 
  matchCompounds <- function(imp_peaks_df, db){
    mass_list_final <- as.numeric(unlist(lapply(imp_peaks_df$mass, 
                                                function(s){strsplit(as.character(s), split = ";")}) ))
    compound_list <- list()
    num_mass_matches <- list()
    big_list <- list()
    MaxZscore <- list()
    i <- 1
    imp_peaks_df_sep <- separate_rows(imp_peaks_df, mass, sep = ";")
    imp_peaks_df_sep$mass <- as.numeric(as.character(imp_peaks_df_sep$mass))
    for (entry in mass_list_final){
      mass_range <- c(entry-(entry*(15/1e6)), entry+(entry*(15/1e6)))
      filtered_db <- db %>% filter(monisotopic_molecular_weight >= mass_range[1], monisotopic_molecular_weight <= mass_range[2])
      compound_list[[i]] <- filtered_db$accession
      num_mass_matches[[i]] <- nrow(filtered_db)
      
      max_zscore <- imp_peaks_df_sep %>% 
        filter(mass == entry)
      MaxZscore[[i]] <- max(as.numeric(max_zscore$zscore), na.rm = TRUE)
      i <- i + 1
    }
    big_list[[1]] <- compound_list
    print("allocated big list 1")
    big_list[[2]] <- num_mass_matches
    print("allocated big list 2")
    big_list[[3]] <- MaxZscore
    print("allocated big list 3")
    return(big_list)
  }
  
  linear_compounds_all_peaks <- matchCompounds(linear_sig_matrix_all_peaks, HMDB)
  print("allocated linear_compounds_all_peaks")
  #background_set_all_peaks are the metabolites that all features (not just sig features) map to
  background_set_all_peaks <- matchCompounds(unn_data, HMDB)
  print("allocated background_compounds_all_peaks")
  
  write.csv(unlist(background_set_all_peaks[[2]]), "number_of_metabolites_mapped_to_features.csv")
  rm(linear_sig_matrix_all_peaks)
  rm(unn_data)
  #prepare data frame with two columns: metabolite HMDB ID and max Z-score found annotated to this metabolite
  prepareCompoundLists <- function(matchCompounds_result, norm_type){
    compound_names_new <- unlist(matchCompounds_result[[1]])
    lengths <- unlist(lapply(matchCompounds_result[[1]], 
                                        function(s){length(s)}) )
    zscores <- as.numeric(unlist(matchCompounds_result[[3]]))
    zscores_all <- rep(zscores, times = lengths)
    names(compound_names_new) <- NULL
    if(norm_type=="background_all"){
      diff_met <- data.frame(met = as.character(compound_names_new))
    }else{
      diff_met <- data.frame(met = compound_names_new, MaxZscore = zscores_all)
    }
    print("allocated diff_met vector")
    diff_met
  }
  
  linear_all_peaks <- prepareCompoundLists(linear_compounds_all_peaks, "linear_all")
  print("allocated linear_all_peaks")
  rm(linear_compounds_all_peaks)
  rm(intensity_matrix)
  rm(peakTable)
  print(ls())
  background_all_peaks <- prepareCompoundLists(background_set_all_peaks, "background_all")
  print("got to right before mode allocation")
 
  rm(background_set_all_peaks)
  if (mode=="positive"){
    linear_all_peaks_pos <- linear_all_peaks
    background_all_peaks_pos <- background_all_peaks
  } else{
    linear_all_peaks_neg <- linear_all_peaks
    background_all_peaks_neg <- background_all_peaks
  }
  
  print("end preparing metabolite lists")
  
  #ENRICHMENT TESTS ---------------------------------
  
  #function for turning csv files into lists
  csv2list <-
  	function (filepath) {
  		
  		# read csv file
  		x <- read.csv(filepath, stringsAsFactors = FALSE)
  		
  		metabolite_set <- as.character(unique(x[,1]))
  		
  		# IDs for metabolite set
  		id <- x[,2]
  		#ID <- NaN
  		ID <- lapply(as.list(metabolite_set), function(set){
  			id_set <- id[x[,1]==set]
  			as.character(unique(id_set))
  		})
  		names(ID) <- as.character(metabolite_set)
  		return(ID)
  	}
  
  #read in list of TIDEX variants from TRIO analysis
  variants <- read.csv(paste0(c("../../candidate_genes/", sample_name, "_variants.csv"), collapse = ""))
  
  #get gene-metabolite annotations

  metabolite_set <- csv2list("../../HMDB_ver4_gene_metabolite_annotations.csv")
  
  MSEA <- function(sig, background, metabolite_set, name_of_file){
    print("in msea")
    names <- names(metabolite_set)
    sig_new <- sig$met
    background_new <- background$met
    pValue <- lapply(metabolite_set, function(x){
      cat <- x
      diff_and_in_cat <- sig_new[sig_new %in% cat]
      diff_and_not_cat <- sig_new[!c(sig_new %in% cat)]
      not_diff_in_cat <- background_new[!c(background_new %in% sig_new)][(background_new[!c(background_new %in% sig_new)]) %in% cat]
      not_diff_not_cat <- background_new[!c(background_new %in% sig_new)][!(background_new[!c(background_new %in% sig_new)]) %in% cat]
      table <- matrix(c(length(diff_and_in_cat), length(diff_and_not_cat), length(not_diff_in_cat),
                        length(not_diff_not_cat)), nrow=2, ncol=2)
      fisher.result <- fisher.test(table, alternative="g")
      fisher.result$p.value
    })
    final_pvalue <- unlist(pValue)
    names(final_pvalue) <- NULL
    qvalue <- p.adjust(final_pvalue, method = "fdr")
    
    MaxZscore <- lapply(metabolite_set, function(x){
      
      diff_and_in_cat <- sig[c(sig_new %in% x),]
      print(diff_and_in_cat)
      print(max(abs(diff_and_in_cat$MaxZscore)))
      if(length(diff_and_in_cat$met) > 0){
        max(abs(diff_and_in_cat$MaxZscore))
      } else{
        0
      }
      
    })
    
    DiffMet <- lapply(metabolite_set, function(x){
      
      paste0(sig[c(sig_new %in% x),]$met, collapse = ", ")
      
    })
    prop <- lapply(metabolite_set, function(x){
      overlap <- intersect(sig_new, x)
      length(overlap)/length(x)
    })
    length <- lapply(metabolite_set, function(x){length(x)})
    
    final <- data.frame(genes = names, 
                        p.value = final_pvalue, 
                        FDR = qvalue, 
                        perc_enrich = unlist(prop), 
                        size = unlist(length), 
                        MaxZscore = unlist(MaxZscore),
                        DAM = unlist(DiffMet)
    )
    write.csv(final, file = name_of_file)
  }

  
  MSEA(linear_all_peaks, background_all_peaks, metabolite_set, paste0(c(mode, "_linear_msea_genebig_allpeaks.csv"), collapse = ""))
  
  print("end msea for single mode")
  setwd("../../")
  getwd()
}


################################################################################################

#combine positive and negative modes and perform MSEA
setwd("./output/")
print("begin combined MSEA")

linear_all_peaks_combined <- rbind(linear_all_peaks_pos, linear_all_peaks_neg)
background_all_peaks_combined <- rbind(background_all_peaks_pos, background_all_peaks_neg)

MSEA(linear_all_peaks_combined, background_all_peaks_combined, metabolite_set, "linear_msea_genebig_allpeaks_combined.csv")

unlink(".RData")
print("end msea")
