library(dplyr)
library(seqinr)
library(tidyr)
library(reshape)
setwd('c:/Users/Alvaro/Documents/Oropouche/Corrected_predictions/')

NetMHCII_parse <- function(data_file, fasta_file, score_type) {
  
  #Generar dataframe
  data <- data.frame(read_delim(data_file, delim = '\t', show_col_types = FALSE, 
                                escape_double = FALSE))
  
  #Obtener longitud y nombre de secuencias
  seqs <- seqinr::read.fasta(fasta_file)
  lenghts <- nchar(getSequence(seqs, as.string = TRUE))-8
  tabla <- data.frame(ID = getName(seqs))
  tabla$length <- lenghts
  data$length <- tabla$length[match(data[, 3], tabla$ID)]
  
  #Obtener nombre de los alelos HLA
  headers <- colnames(data)
  HLA_alleles <- headers[grep("DRB", headers)]
  
  #Indicar qué columnas se extraen según el score_type
  if (score_type == "Rank") {
    cols_extract <- seq(from = 7, by = 3, length.out = length(HLA_alleles))
  } else {
    cols_extract <- seq(from = 6, by = 3, length.out = length(HLA_alleles))
  }
  
  #Crear el archivo de salida
  NetMHCII_parsed <- data[-1, c(2, 1, ncol(data), 3, cols_extract)]
  colnames(NetMHCII_parsed) <- c("Peptide", "Position", "Length", "ID", 
                                 HLA_alleles)
  NetMHCII_parsed$Position <- as.numeric(NetMHCII_parsed$Position)
  NetMHCII_parsed[c(5:ncol(NetMHCII_parsed))] <- 
    sapply(NetMHCII_parsed[c(5:ncol(NetMHCII_parsed))], as.numeric)
  NetMHCII_parsed$Position <- as.numeric(NetMHCII_parsed$Position)

  return(NetMHCII_parsed)
}

MixMHC2pred_parse <- function(prediction_file, auxiliary_file, fasta_file){
  
  #Generar dataframe
  data <- data.frame(read_delim(prediction_file, delim = "\t", 
                                   show_col_types = FALSE, escape_double = FALSE))
  
  #Obtener una tabal de donde se van a sacar las longitudes, ID y posiciones
  auxiliary <- function(auxiliary_file, fasta_file) {
    data_auxiliary <- data.frame(read_delim(auxiliary_file, delim = '\t', 
                                            show_col_types = FALSE, 
                                            escape_double = FALSE))
    seqs <- seqinr::read.fasta(fasta_file)
    lenghts <- nchar(getSequence(seqs, as.string = TRUE))-8
    tabla <- data.frame(ID = getName(seqs))
    tabla$length <- lenghts
    data_auxiliary$length <- tabla$length[match(data_auxiliary[, 3], tabla$ID)]
    return(data_auxiliary)
  }
  
  data_auxiliary <- auxiliary(auxiliary_file, fasta_file)
  
  #Añadir ID, Position y length a la data
  data$ID <- data_auxiliary$...3[match(data[, 1], 
                                             data_auxiliary$...2)]
  data$length <- data_auxiliary$length[match(data[, 1], 
                                             data_auxiliary$...2)]
  data$Position <- data_auxiliary$...1[match(data[, 1], 
                                             data_auxiliary$...2)]
  
  #Extraer columnas con el Rank
  HLA_alleles <- data[grep("Rank_DRB", colnames(data))]
  
  #Generar tabla con resultados
  pept_features <- data[, c(1, ncol(data), ncol(data)-2, ncol(data)-1)]
  MixMHC2_parse <- bind_cols(pept_features, HLA_alleles, id = NULL)
  MixMHC2_parse$Position <- as.numeric(MixMHC2_parse$Position)
  MixMHC2_parse[c(5, ncol(MixMHC2_parse))] <- 
    sapply(MixMHC2_parse[c(5, ncol(MixMHC2_parse))], as.numeric)
  
  return(MixMHC2_parse)
}


Merge_results_MHC2 <- function(NetMHC2_file, MixMHC2_file, fasta_file, score_type) {
  
  #Crear archivos ordenados
  file1 = NetMHCII_parse(NetMHC2_file, fasta_file, score_type)
  file2 = MixMHC2pred_parse(MixMHC2_file, NetMHC2_file, fasta_file)
  
  file1 = file1 %>% arrange(ID, Position)
  file2 = file2 %>% arrange(ID, Position)
  
  #Crear output mergeado
  Merged_SudAm_MHC2 = file1
  
  for (i in 1:nrow(Merged_SudAm_MHC2)) {
    for (j in 5:ncol(Merged_SudAm_MHC2)) {
      if (file1[i, j] > 5 | file2[i, j] > 5) {
        Merged_SudAm_MHC2[i, j] = NA
      }
    }
  }
  
  Merged_SudAm_MHC2$Position = as.character(Merged_SudAm_MHC2$Position)
  Merged_SudAm_MHC2$Length = as.character(Merged_SudAm_MHC2$Length)
  
  Merged_SudAm_MHC2 = Merged_SudAm_MHC2 %>%
    filter(!rowSums(is.na(select_if(., is.numeric))) == ncol(select_if(., is.numeric)))
  
  write.table(Merged_SudAm_MHC2, 'Procesed/Merged_Global_MHC2.txt', sep = '\t', 
              row.names = FALSE, quote = FALSE)
  
}

Merge_results_MHC2('New_Results/NetMHCIIpan_Results_Global.txt', 'New_Results/MixMHC2_Results_Global.txt',
                   'Proteins_consenso.fasta', 'Rank')


file1 = NetMHCII_parse('NetMHCIIpan/NetMHCIIpan_Global.xls', 'Gn_NSm_Proteins_consenso.fasta', 'Rank')
file2 = MixMHC2pred_parse('MixMHC2pred/Gn_NSm_MixMHC_Global.txt.txt', 'NetMHCIIpan/NetMHCIIpan_Global.xls', 
                          'Gn_NSm_Proteins_consenso.fasta')

file1 = file1 %>% arrange(ID, Position)
file2 = file2 %>% arrange(ID, Position)

file2[c(5, ncol(file2))] = 
  sapply(file2[c(5, ncol(file2))], as.numeric)

Merged_Global_MHC2 = file1
for (i in 1:nrow(Merged_Global_MHC2)) {
  for (j in 5:ncol(Merged_Global_MHC2)) {
    if (file1[i,j] > 5 | file2[i,j] > 5) {
      Merged_Global_MHC2[i,j] = NA
    }
  }
}

Merged_Global_MHC2$Position = as.character(Merged_Global_MHC2$Position)
Merged_Global_MHC2$Length = as.character(Merged_Global_MHC2$Length)

Merged_Global_MHC2 = Merged_Global_MHC2 %>%
  filter(!rowSums(is.na(select_if(., is.numeric))) == ncol(select_if(., is.numeric)))

write.table(Merged_Global_MHC2, 'Processed/Gn_NSm_Merged_Global_MHC2.txt', sep = '\t', 
            row.names = FALSE, quote = FALSE)

#Converitr NA a 9
Merged_Global_MHC2[is.na(Merged_Global_MHC2)] <- 9

write.table(Merged_Global_MHC2, 'Processed_with9/Merged_SudAm_MHC2_na9.txt', 
            sep = "\t",row.names = FALSE, quote = FALSE)
