library(dplyr)
library(readr)
library(reshape)
library(seqinr)
library(tidyr)
setwd('c:/Users/Alvaro/Documents/Oropouche/Corrected_predictions/')

Parse_NetMHCpan <- function(prediction_file, fasta_file, score_type) {
  
  #Crear dataframe
  data <- data.frame(read_delim(prediction_file, delim = '\t', 
                               escape_double = FALSE, show_col_types = FALSE))
  
  #Obtener longitudes de secuencias
  fasta <- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(fasta, as.string = TRUE))-8
  
  #Obtener nombres de secuencias
  seq.names <- data.frame(names = getName(fasta))
  table <- tidyr::separate(seq.names, names, into ="ID")
  
  #Añadir longitudes de las secuencia a esta tabla
  table$length <- lengths
  
  #Obtener los alelos HLA
  headers <- colnames(data)
  alleles_hla <- headers[grep("HLA", headers)]
  
  #Añadir longitudes a la data
  data$length <- table$length[match(data[,3], table$ID)]
  
  #Indicar de qué columnas se va a extraer el score de cada epítopo para cada alelo
  if (score_type == "Rank") {
    col_extracted <- seq(from = 7, by = 4, length.out = length(alleles_hla))
  } else {
    col_extracted <- seq(from = 6, by = 4, length.out = length(alleles_hla))
  }
  
  #Generar el dataframe output
  NetMHC_parsed <- data[-1, c(2, 1, ncol(data), 3, col_extracted)]
  colnames(NetMHC_parsed) <- c("Peptide", "Position", "Length", "ID", alleles_hla)
  NetMHC_parsed$Position <- as.numeric(NetMHC_parsed$Position)+1
  NetMHC_parsed[c(5:ncol(NetMHC_parsed))] <- 
    sapply(NetMHC_parsed[c(5:ncol(NetMHC_parsed))], as.numeric)
  
  return(NetMHC_parsed)
}


Parse_MHCFlurry <- function(prediction_file, fasta_file, score_type) {
  
  #Crear dataframe
  data <- data.frame(read_delim(prediction_file, delim = ',', 
                                escape_double = FALSE, show_col_types = FALSE))
  
  #Obtener longitud de secuencias
  fasta <- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(fasta, as.string = TRUE))-8
  
  #Obtener nombres de secuencias
  seq.names <- data.frame(names = getName(fasta))
  table <- tidyr::separate(seq.names, names, into = "ID")
  
  #Añadir longitud de secuencias
  table$length <- lengths
  
  #Ordenar al formato adecuado los resultados de MHCFlurry
  if (score_type == "Rank") {
    data <- data[, c(1, 2, 3, 8, 9)]
  } else {
    data <- data[, c(1, 2, 3, 8, 7)]
  }
  arrange <- cast(data, formula = sequence_name + pos + peptide ~ best_allele)
  
  #Añadir longitud a la tabla
  arrange$length <- table$length[match(arrange[, 1], table$ID)]
  
  #Generar el dataframe output
  MHCFlurry_parsed <- arrange[, c(3, 2, ncol(arrange), 1, 4:(ncol(arrange)-1))]
  colnames(MHCFlurry_parsed)[1:4] <- c("Peptide", "Position", "Length", "ID")
  MHCFlurry_parsed$Position <- as.numeric(MHCFlurry_parsed$Position)+1
  MHCFlurry_parsed[c(5:ncol(MHCFlurry_parsed))] <- 
    sapply(MHCFlurry_parsed[c(5:ncol(MHCFlurry_parsed))], as.numeric)
  
  return(MHCFlurry_parsed)
}

#Ejecuta con los archivos aquí abajo umu
file1 = Parse_NetMHCpan('NetMHCpan/NetMHCpan_Global.xls', 'Gn_NSm_Proteins_consenso.fasta',
                        'Rank')
file2 = Parse_MHCFlurry('MHCFlurry/Global_epitopes_mhcflurry.csv', 'Gn_NSm_Proteins_consenso.fasta',
                        'Affinity')

file1 = file1 %>% arrange(ID, Position, Peptide)
file2 = file2 %>% arrange(ID, Position, Peptide)

#Unión de los resultados de NetMHC y MHCFlurry
Merged_Global_MHC1 = file1
for (i in 1:nrow(Merged_Global_MHC1)) {
  for (j in 5:ncol(Merged_Global_MHC1)) {
    if (file1[i, j] > 2 | file2[i, j] > 500) {
      Merged_Global_MHC1[i, j] = NA
    }
  }
}

Merged_Global_MHC1$Position = as.character(Merged_Global_MHC1$Position)
Merged_Global_MHC1$Length = as.character(Merged_Global_MHC1$Length)


Merged_Global_MHC1 = Merged_Global_MHC1 %>%
  filter(!rowSums(is.na(select_if(., is.numeric))) == ncol(select_if(., is.numeric)))

write.table(Merged_Global_MHC1, 'Processed/Gn_NSm_Merged_Global_MHC1.txt', sep = '\t', 
            row.names = FALSE, quote = FALSE)

#Converitr NA a 9
Merged_Global_MHC1[is.na(Merged_Global_MHC1)] <- 9

write.table(Merged_Global_MHC1, 'Processed_with9/Merged_SudAm_MHC1_na9.txt', 
            sep = "\t",row.names = FALSE, quote = FALSE)




