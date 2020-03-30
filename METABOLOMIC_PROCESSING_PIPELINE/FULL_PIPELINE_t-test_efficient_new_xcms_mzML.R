#FULL PIPELINE

#LOAD LIBRARIES
unlink(".RData")
options(BIOCINSTALLER_ONLINE_DCF = FALSE)
.libPaths("/zfs3/users/grahamemma9/grahamemma9/R/3.5/lib64/R/library")
library(norm)
library(xcms, lib.loc = "/zfs3/users/grahamemma9/grahamemma9/R/x86_64-pc-linux-gnu-library/3.5")
library(CAMERA, lib.loc = "/zfs3/users/grahamemma9/grahamemma9/R/3.5/lib64/R/library")
library(XML, lib.loc = "/zfs3/users/grahamemma9/grahamemma9/R/3.5/lib64/R/library")
library(reshape2, lib.loc = "/zfs3/users/grahamemma9/grahamemma9/R/3.5/lib64/R/library")
library(tidyverse, lib.loc = "/zfs3/users/grahamemma9/grahamemma9/R/3.5/lib64/R/library")
#library(mseapca, lib.loc = "/zfs3/users/grahamemma9/grahamemma9/R/3.5/lib64/R/library")
args=(commandArgs(TRUE))
sample_name <-as.character(args[[1]])

#XCMS------------------------------------


for (mode in c("positive", "negative")){
    
  if(mode == "positive"){
    control_files <- list.files(path = "Controls/Files_pos", recursive = TRUE, full.names = TRUE)
    patient_files <- list.files(path = paste0(c(sample_name, "/Files_pos"), collapse = ""), recursive = TRUE, full.names = TRUE)
  } else{
    control_files <- list.files(path = "Controls/Files_neg", recursive = TRUE, full.names = TRUE)
    patient_files <- list.files(path = paste0(c(sample_name, "/Files_neg"), collapse = ""), recursive = TRUE, full.names = TRUE)
  }
  print(control_files)
  print(patient_files)
  
  our_files <-  c(control_files, patient_files)
  
  #pick peaks
  meta_data <- data.frame(sample_name = sub(basename(our_files), pattern = ".mzML",
  																					replacement = "", fixed = TRUE),
  												sample_group = c(rep("Control", 10), sample_name),
  												stringsAsFactors = FALSE) 
  if(file.exists(paste0(c("./",sample_name, "/",sample_name,"_", mode, "_xsaFA.RDS"), collapse = "")) == TRUE){
    print("yes")
    xsaFA <- readRDS(paste0(c("./",sample_name, "/",sample_name,"_", mode, "_xsaFA.RDS"), collapse = ""))
    print("got past load RDS")
  }else{
    print("no")
    raw_data <- readMSData(files = our_files, pdata = new("NAnnotatedDataFrame", meta_data),
  											 mode = "onDisk") 
  # xset <- xcmsSet(our_files,
  # 								method = "centWave",
  # 								ppm = 15,
  # 								peakwidth = c(3,35.75),
  # 								mzdiff = 0.00325,
  # 								prefilter = c(3, 100),
  # 								noise = 0,
  # 								snthresh = 2.8)
  #cwp <- CentWaveParam(peakwidth = c(30, 80), noise = 1000)
    cwp <- CentWaveParam(peakwidth = c(3, 80), 
  										 noise = 1000,
  										 ppm = 15,
  										 snthresh = 10,
  										 mzdiff = 0.0045,
  										 prefilter = c(3, 100))
  xdata <- findChromPeaks(raw_data, param = cwp)
  pdp <- PeakDensityParam(sampleGroups = xdata$sample_group)
  xdata <- groupChromPeaks(xdata, param = pdp)
  xdata <- adjustRtime(xdata, param = ObiwarpParam(gapInit = 1.2736,
  																								 gapExtend = 3.3336))
  pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
  												minFraction = 1, 
  												bw = 0.25,
  												minSamples = 1)
  xdata <- groupChromPeaks(xdata, param = pdp)
  any(is.na(chromPeaks(xdata)))
  xdata <- fillChromPeaks(xdata)
  
  intensity_matrix <- featureValues(xdata, method = "medret", value = "into")
  feature_def <- featureDefinitions(xdata)
  write.csv(intensity_matrix, paste0(c("./", sample_name, "/", mode, "_intensity_matrix_April2019_medret_into.csv"), collapse = ""))
  write.csv(feature_def, paste0(c("./", sample_name, "/", mode, "_feature_definitions_April2019_medret_into.csv"), collapse = ""))

  print("got to xcmsSet")
  
  #CAMERA-------------------------------------
  
  #setwd("/Users/emmagraham/Desktop/Masters/Metabolomics project/singlePatient_vs_allCTRLs/TX071_patient metabolome")
  #Create an xsAnnotate object
  xcmsSetObject <- as(xdata, "xcmsSet")
  xcmsSetObject <- fillPeaks(xcmsSetObject)
  xsa <- xsAnnotate(xcmsSetObject, sample = c(1:length(our_files)), polarity = mode)
  #Group after RT value of the xcms grouped peak
  xsaF <- groupFWHM(xsa, perfwhm=0.6)
  #Verify grouping
  xsaC <- groupCorr(xsaF)
  #Annotate isotopes, could be done before groupCorr
  xsaFI <- findIsotopes(xsaC)
  #Annotate adducts
  xsaFA <- findAdducts(xsaFI, polarity=mode)
  saveRDS(xsaFA, file = paste0(c("./",sample_name, "/",sample_name,"_", mode, "_xsaFA.RDS"), collapse = ""))
  }
  library(tidyverse)
  #Get final peaktable and store
  print("xsaFA")
  peaks <- getPeaklist(xsaFA, intval = "intb") %>% tbl_df()
  
  colnames(peaks) <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", 
                                                       "Control", sample_name, unlist(lapply(control_files, function(x){substring(x, nchar(x)-16, nchar(x))})), unlist(lapply(patient_files, function(x){substring(x, nchar(x)-16, nchar(x))})), 
                                                      "isotopes",
                                                      "adduct", "pcgroup")
  peakTable_columns <- peaks %>% 
    dplyr::select(mz, mzmin, mzmax, rt, rtmin, rtmax, npeaks, Control, sample_name, isotopes, adduct, pcgroup)
  intensity_matrix <- read.csv(paste0(c("./", sample_name, "/", mode, "_intensity_matrix_April2019_medret_into.csv"), collapse = ""),
                               stringsAsFactors = FALSE, 
                               row.names = 1)
  print(c(unlist(lapply(control_files, function(x){substring(x, nchar(x)-16, nchar(x))})), unlist(lapply(patient_files, function(x){substring(x, nchar(x)-16, nchar(x))}))))
  peakTable <- cbind(peakTable_columns, intensity_matrix[, c(unlist(lapply(control_files, function(x){substring(x, nchar(x)-16, nchar(x))})), unlist(lapply(patient_files, function(x){substring(x, nchar(x)-16, nchar(x))})))])
  dim(peakTable)
  peakTable <- peakTable[is.finite(rowSums(peakTable[, c(unlist(lapply(control_files, function(x){substring(x, nchar(x)-16, nchar(x))})),
                                       unlist(lapply(patient_files, function(x){substring(x, nchar(x)-16, nchar(x))})))])),]
  dim(peakTable)
  write.csv(peakTable, paste0(c("./", sample_name, "/", mode, "peakTable.csv"), collapse = ""))
  print("done CAMERA")
  #STATS
  
  setwd(paste0(c("./", sample_name), collapse = ""))
  if (mode == "positive"){
    if(dir.exists("stuff_pos")){
      setwd("./stuff_pos")
    } else{
      dir.create("stuff_pos")
      setwd("./stuff_pos")
    }
  } else{
    if(dir.exists("stuff_neg")){
      setwd("./stuff_neg")
    } else{
      dir.create("stuff_neg")
      setwd("./stuff_neg")
    }
  }
  
  #get data frame of intensities from DIFFREPORT
  intensity_matrix <- read.csv(paste0(c("../", mode, "_intensity_matrix_April2019_medret_into.csv"), collapse = ""),
                                stringsAsFactors = FALSE, 
                               row.names = 1)
  write.csv(intensity_matrix, "./unnormalized_raw_intensities_April_2019.csv")
  unnormalized_data <- intensity_matrix[is.finite(rowSums(log2(intensity_matrix))),]
  log_data <- log2(unnormalized_data)[is.finite(rowSums(log2(unnormalized_data))),]
  write.csv(log_data, "./log_raw_intensities.csv")
  
  #linear baseline
  
  linear.baseline <- apply(log_data,1, function(x){median(x[complete.cases(x)])}) 
  print(linear.baseline)
  baseline.mean <- mean(linear.baseline, na.rm = TRUE)
  print(baseline.mean)
  sample.means <- colMeans(log_data, na.rm = TRUE) 
  print(sample.means)
  linear.factor <- baseline.mean/sample.means 
  print(linear.factor)
  linear.normalized.data <- t(t(log_data)*linear.factor)
  write.csv(linear.normalized.data, "./linear_raw_intensities_April_2019.csv")
  
  print("done with normalization")
  
  #remove lowly expressed  peaks i.e those that are in less than half of control and 
  #less than half of patient. These are likely not important peaks. 
  library(tidyverse)
  # identifyBadPeaks <- function(data){
  #   discard <- c()
  #   for (entry in 1:nrow(data)){
  #     if (as.data.frame(data)[, 8][entry] <= 7 && as.data.frame(data)[, 9][entry] == 0){
  #       discard <- c(discard, entry)
  #     } else{
  #       0
  #     }
  #   }
  #   return(discard)
  # }
  # badPeaksIndex <- identifyBadPeaks(peakTable)
  # print("done with finding bad peaks")
  # #removes bad peaks
  # applySelection <- function(data, index){
  #   data_wo_bad_peaks <- data[-c(index),]
  #   return(data_wo_bad_peaks)
  # }
  
  #this function will generate a table with a mz, rt, isotope and adduct column. The input is the 
  #intensity matrix,  an index of peaks to remove and the table with adduct and isotope annotation.
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
  print("linear_normalized_w_columns")
  linear_normalized_w_columns <- addColumns(linear.normalized.data, badPeaksIndex, peakTable, FALSE)
  str(linear_normalized_w_columns)
  unnormalized_w_columns <- addColumns(unnormalized_data, badPeaksIndex, peakTable, FALSE)
  
  
  #Isolate annotated MM from adduct column. If a compound could be multiple MM, all are included. 
  #If no adduct is annotated, then the unannotated peak mz is used in its place. 
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
  #Add masses to data wo low peaks
  print("linear_normalized_w_columns_w_masses")
  linear_normalized_w_columns_w_masses <- addMasses(linear_normalized_w_columns)
  print(str(linear_normalized_w_columns_w_masses))
  unnormalized_w_columns_w_masses <- addMasses(unnormalized_w_columns)
  
  #this function removes isotopes that are not base isotopes = [M]+
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
  
  print("linear_normalized_w_columns_w_masses_wo_iso")
  linear_normalized_w_columns_w_masses_wo_iso <- removeIsotopes(linear_normalized_w_columns_w_masses)
  print(str(linear_normalized_w_columns_w_masses_wo_iso))
  unnormalized_w_columns_w_masses_wo_iso <- removeIsotopes(unnormalized_w_columns_w_masses)
  
  print("done data cleaning")
  
  #ISOLATE FEATURES THAT MAP TO COMPOUNDS
  #read in csv file containing HMDB database
  
  HMDB <- read.csv("../../HMDB_db_ver4.csv", header = TRUE)
  print(str(HMDB))
  
  #isolate control and patient names
  #setwd("../..")
  control_names <- unlist(lapply(control_files, function(x){substring(x, nchar(x)-16, nchar(x))}))
  patient_names <- unlist(lapply(patient_files, function(x){substring(x, nchar(x)-16, nchar(x))}))
  print(c(control_names, patient_names))
  #setwd(paste0(c("./", sample_name, "/wilcox_ebam_fixed_p0"), collapse = ""))
  
  linear_data <- linear_normalized_w_columns_w_masses_wo_iso
  unn_data <- unnormalized_w_columns_w_masses_wo_iso
  
  #function that matches metabolite names with masses within 15ppm. 
  
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
  
  print("linear_compounds_all_peaks1")
  linear_compounds_all_peaks1 <- matchCompounds(linear_data, HMDB)
  print(str(linear_compounds_all_peaks1))
  unn_compounds_all_peaks1 <- matchCompounds(unn_data, HMDB)
  
  print("end preparing metabolite lists")
  
  #WILCOXON-RANK SUM TEST or other statistical tests----------------
  
  
  linear_data1 <- t(linear_compounds_all_peaks1[,colnames(linear_compounds_all_peaks1) %in% c(control_names, patient_names)])
  scaled <- apply(linear_data1, 2, function(x){scale(x, center = TRUE, scale = TRUE)})
  special <- abs(scaled[11, ]) >= 2
  print(dim(linear_compounds_all_peaks1))
  linear_compounds_all_peaks1$zscore <- scaled[11,]
  linear_sig_matrix_all_peaks <- linear_compounds_all_peaks1[special,]
  write.csv(linear_sig_matrix_all_peaks, "./significant_features_linear_April_2019.csv")
  print("done with stats")
  #HMDB lookup-------------------------------
  
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
      #ADD HERE 
      i <- i + 1
    }
    big_list[[1]] <- compound_list
    big_list[[2]] <- num_mass_matches
    big_list[[3]] <- MaxZscore
    return(big_list)
  }
  
  linear_compounds_all_peaks <- matchCompounds(linear_sig_matrix_all_peaks, HMDB)
  background_set_all_peaks <- matchCompounds(unn_data, HMDB)
  print(memory.size(background_set_all_peaks))
  write.csv(unlist(background_set_all_peaks[[2]]), "features_mapped_to_compounds_April_2019.csv")
  
  prepareCompoundLists <- function(matchCompounds_result, norm_type){
    compound_names_new <- as.character(unlist(matchCompounds_result[[1]]))
    print(compound_names_new)
    lengths <- unlist(lapply(matchCompounds_result[[1]], 
                             function(s){length(s)}) )
    zscores <- as.numeric(unlist(matchCompounds_result[[3]]))
    zscores_all <- rep(zscores, times = lengths)
    names(compound_names_new) <- NULL
    diff_met <- data.frame(met = compound_names_new, MaxZscore = zscores_all)
    write.csv(diff_met, file = "sig_compound_names.csv",
              quote = FALSE,
              row.names = FALSE)
    print("got to right before diff met output")
    diff_met
  }
  linear_all_peaks <- prepareCompoundLists(linear_compounds_all_peaks, "linear_all")
  #background_all_peaks <- prepareCompoundLists(background_set_all_peaks, "background_all")
  print("allocated linear all peaks vector")
  print(length(as.character(unlist(background_set_all_peaks[[1]]))))
  background_all_peaks <- as.character(unlist(background_set_all_peaks[[1]]))
  print("got to background_all_peaks")
  if (mode=="positive"){
    linear_all_peaks_pos <- linear_all_peaks
    background_all_peaks_pos <- background_all_peaks
  } else{
    linear_all_peaks_neg <- linear_all_peaks
    background_all_peaks_neg <- background_all_peaks
  }
  
  print("end preparing metabolite lists")
  
  #get variant list
  
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
      # for (i in 1:length(metabolite_set)){
      #   id_set <- id[x[,1]==metabolite_set[i]]
      #   ID[i] <- list(unique(id_set))
      # }
      # 
      # names(ID) <- metabolite_set
      # return(ID)
    }
  
  
  variants <- read.csv(paste0(c("../../TIDEX_variants/", sample_name, "_variants.csv"), collapse = ""))
  
  #get gene-metabolite annotations
  metabolite_set_list <- read.csv("../../full_metabolite_set_HMDB_ver4_outofthebox.csv")
  final_metabolite_set <- metabolite_set_list %>% 
    filter(gene_name %in% variants$gene_name)
  
  #metabolite_set_singleton <- read.csv("../../full_metabolite_gene_set.csv")
  #final_metabolite_set <- metabolite_set %>% 
  #	filter(metabolite.set.name %in% variants_singleton$gene_name)
  
  # final_metabolite_set <- data.frame(metabolite.set.name = NA, metabolite.IDs = NA)
  # for (entry in 1:nrow(metabolite_set_small)){
  #   name <- metabolite_set_small$metabolite.set.name[entry]
  #   compounds <- unlist(lapply(metabolite_set_small$metabolite.IDs[entry],
  #                              function(s){strsplit(as.character(s), split = "; ")}) )
  #   if (length(compounds)>0){
  #     foo <- lapply(compounds, function(x){
  #       data.frame(metabolite.set.name = name, metabolite.IDs = x)
  #     })
  #     one_set <- do.call("rbind", foo)
  #     final_metabolite_set <- rbind(final_metabolite_set, one_set)
  #   } else{
  #     0
  #   }
  #   print(entry)
  # }
  
  write.csv(final_metabolite_set[complete.cases(final_metabolite_set),], file = "../metabolite_set_small_final.csv", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  metabolite_set_list_small <- csv2list("../metabolite_set_small_final.csv")
  metabolite_set <- csv2list("../../full_metabolite_set_HMDB_ver4_outofthebox.csv")
  MSEA <- function(sig, background, metabolite_set, name_of_file){
    print("in msea")
    names <- names(metabolite_set)
    sig_new <- sig$met
    background_new <- background
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
  
  #in_SMPDB <- read.csv("../../metabolites 5.csv")$HMDB.ID
  
  # MSEA_pathway <- function(sig, background, metabolite_set, name_of_file){
  #   names <- names(metabolite_set)
  #   sig_new <- sig[sig %in% in_SMPDB]
  #   background_new <- background[background %in% in_SMPDB]
  #   pValue <- lapply(metabolite_set, function(x){
  #     cat <- x
  #     diff_and_in_cat <- sig_new[sig_new %in% cat]
  #     diff_and_not_cat <- sig_new[!c(sig_new %in% cat)]
  #     not_diff_in_cat <- background_new[!c(background_new %in% sig_new)][(background_new[!c(background_new %in% sig_new)]) %in% cat]
  #     not_diff_not_cat <- background_new[!c(background_new %in% sig_new)][!(background_new[!c(background_new %in% sig_new)]) %in% cat]
  #     table <- matrix(c(length(diff_and_in_cat), length(diff_and_not_cat), length(not_diff_in_cat),
  #                       length(not_diff_not_cat)), nrow=2, ncol=2)
  #     fisher.result <- fisher.test(table, alternative="g")
  #     fisher.result$p.value
  #   })
  #   final_pvalue <- unlist(pValue)
  #   names(final_pvalue) <- NULL
  #   qvalue <- p.adjust(final_pvalue, method = "fdr")
  #   
  #   prop <- lapply(metabolite_set, function(x){
  #     overlap <- intersect(sig_new, x)
  #     length(overlap)/length(x)
  #   })
  #   length <- lapply(metabolite_set, function(x){length(x)})
  #   final <- data.frame(genes = names, p.value = final_pvalue, FDR = qvalue, perc_enrich = unlist(prop), size = unlist(length))
  #   write.csv(final, file = name_of_file)
  # }
  
  #gene-metabolite sets
  #all peaks gene metabolite lists
  
  MSEA(linear_all_peaks, background_all_peaks, metabolite_set_list_small, paste0(c(mode, "_linear_allpeaks_msea_genesmall.csv"), collapse = ""))
  
  #metabolism sets
  
  #pathway_set_list <- csv2list("../../pathway_sets.csv")
  
  
  #MSEA_pathway(linear_all_peaks, background_all_peaks, pathway_set_list, paste0(c(mode, "_linear_msea_pathway_all_peaks.csv"), collapse = ""))
  
  #gene-metabolite sets with full gene list
  print(class(metabolite_set))
  
  MSEA(linear_all_peaks, background_all_peaks, metabolite_set, paste0(c(mode, "_linear_msea_genebig_allpeaks_April_2019.csv"), collapse = ""))
  
  #singleton variants
  if (sample_name %in% c("VN0149", "TX400", "TIDEX748", "VN0121", "TIDEX698", "TIDEX540",
                         "VN0099", "VB0028", "TX071", "TX039", "TIDEX042", "TX196", "TIDEX735")){
    variants_singleton <- read.csv(paste0(c("../../TIDEX_variants/Oct31_update/", sample_name, "_Oct31_update.csv"), collapse = ""),
                                   stringsAsFactors = FALSE)
    final_metabolite_set_singleton <- metabolite_set_list %>% 
      filter(gene_name %in% variants_singleton$gene_name)
    
    
    write.csv(final_metabolite_set_singleton[complete.cases(final_metabolite_set_singleton),], 
              file = "../metabolite_set_small_final_singleton.csv", 
              quote = FALSE, 
              row.names = FALSE, 
              col.names = TRUE)
    
    metabolite_set_list_small_singleton <- csv2list("../metabolite_set_small_final_singleton.csv")
    MSEA(linear_all_peaks, background_all_peaks, metabolite_set_list_small_singleton, 
         paste0(c(mode, "_linear_allpeaks_msea_genesmall_singleton_April_2019.csv"), collapse = ""))
    #MSEA(linear_all_peaks, background_all_peaks, metabolite_set_list, paste0(c(mode, "_linear_msea_genebig_allpeaks_singleton.csv"), collapse = ""))
  } else{
    0
  }
  
  print("end msea")
  setwd("../..")
  getwd()
}



# par(mfrow = c(3, 2))
# hist(melt(unnormalized_data[1:1000,])$value, breaks = 1000, xlim = c(0, 1e5))
# hist(melt(linear.normalized.data[1:1000,])$value, breaks = 1000, xlim = c(0, 100))
# hist(melt(loess.data[1:1000,])$value, breaks = 1000, xlim = c(0, 4e4))
# hist(melt(log_data[1:1000,])$value, breaks = 1000, xlim = c(0, 20))
# hist(melt(TIC_normalized[1:1000,])$value, breaks = 1000, xlim = c(0, 4e4))


################################################################################################

#COMBINED MSEA
setwd(paste0(c("./", sample_name), collapse = ""))
linear_all_peaks_combined <- rbind(linear_all_peaks_pos, linear_all_peaks_neg)
background_all_peaks_combined <- rbind(background_all_peaks_pos, background_all_peaks_neg)

MSEA(linear_all_peaks_combined, background_all_peaks_combined, metabolite_set_list_small, "linear_allpeaks_msea_genesmall_combined_April_2019.csv")

#metabolism sets

#MSEA_pathway(linear_all_peaks_combined, background_all_peaks_combined, pathway_set_list, "linear_msea_pathway_all_peaks_combined.csv")

#gene-metabolite sets with full gene list

MSEA(linear_all_peaks_combined, background_all_peaks_combined, metabolite_set, "linear_msea_genebig_allpeaks_combined_April_2019.csv")
if (sample_name %in% c("VN0149", "TX400", "TIDEX748", "VN0121", "TIDEX698", "TIDEX540",
                       "VN0099", "VB0028", "TX071", "TX039", "TIDEX042", "TX196", "TIDEX735")){
  MSEA(linear_all_peaks_combined, background_all_peaks_combined, metabolite_set_list_small_singleton, "linear_allpeaks_msea_genesmall_combined_singleton_April_2019.csv")
} else{0}

unlink(".RData")
print("end msea")
