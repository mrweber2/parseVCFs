### ParseVCFs: functions for file conversion and mutation frequency calculation

## ConvertVCF: function for converting vcf to bed and CSV formats (Default = bed)
##             Returns bed or CSV file
convertVCF <- function(vcf, format = "bed") {
  ## Read vcf lines
  con <- file(vcf, open = 'r')
    
  line_count = 0
    
  ## While loop to read vcf connection
  while(length(lines <- readLines(con))) {
    ## Define columns
    fields <- strsplit(lines, '\t')
    start_line <- 0
    
    ## Count useable lines, skipping the header
    IS_FIRST <- TRUE
    for(i in 1:length(lines)) {
      if(grepl("#", fields[[i]]) == FALSE) {
        if(IS_FIRST == TRUE) {
          start_line <- i  ## Grab line number of first non-header line
          line_count <- line_count + 1
          IS_FIRST <- FALSE
        } else {
          line_count <- line_count + 1
        }
      }
    }
      
    ## Initialize new file, length of line count plus 1 for new header
    new_file <- vector(length = line_count+1)
    
    total_length <- start_line + line_count
    
    ## Add converted lines to new file
    for(i in start_line:length(fields)) {
      ## Grab only the chromosome, start position, alt/ref, and genotype fields
      ## Save to vector at index=line number
      Chrom <- fields[[i]][1]
      Pos <- fields[[i]][2]
      Ref <- fields[[i]][4]
      Alt <- fields[[i]][5]
      GT <- fields[[i]][10]
      
      ## Convert to tab-separated bed file or CSV with Chromosome, position, allele change, and genotype
      if(grepl("bed", format, fixed = TRUE) == TRUE | grepl("BED", format, fixed = TRUE) == TRUE) {
        new_file[(start_line+2)-i] <- paste(Chrom, Pos, Ref, Alt, GT, sep = "\t")
      } else
        if(grepl("csv", format, fixed = TRUE) == TRUE | grepl("CSV", format, fixed = TRUE) == TRUE) {
          new_file[(start_line+2)-i] <- paste(Chrom, Pos, Ref, Alt, GT, sep = ",")
        }
    }
  }
  
  ## Set header as index 1/first line of new file
  if(grepl("bed", format, fixed = TRUE) == TRUE | grepl("BED", format, fixed = TRUE) == TRUE) {
    new_file[1] <- "CHROM\tPOS\tREF\tALT\tGT"
  } else
    if(grepl("csv", format, fixed = TRUE) == TRUE | grepl("CSV", format, fixed = TRUE) == TRUE) {
      print("WARN: Standard VCF contains commas in INFO fields.\nConverted file will have additional commas in the genotype field.")
      new_file[1] <- "CHROM,POS,REF,ALT,GT"
    }
  
  return(new_file)
  
}

## OverallMutRates: function to calculate mutation frequencies for point mutations. 
##                  Returns vector list of mutation frequency with names = allele changes
overallMutRates <- function(vcf) {

  ## Open vcf connection for reading
  con <- file(vcf, open = 'r')

  ## Define possible point mutations
  A2C <- 0
  A2G <- 0
  A2T <- 0
  C2A <- 0
  C2G <- 0
  C2T <- 0
  G2A <- 0
  G2C <- 0
  G2T <- 0
  T2A <- 0
  T2C <- 0
  T2G <- 0
  
  ## While loop to read VCF connection
  while(length(lines <- readLines(con))) {
    ## Skip header
    header <- substr(lines, 1, 1) == "#"
    
    ## Define columns
    fields <- strsplit(lines, '\t')
    
    ## Initialize reference and alt allele vectors
    Ref <- vector(length = length(lines))
    Alt <- vector(length = length(lines))
    
    ## Save Alt and Ref alleles to lists
    for(i in 1:length(lines)) {
      Ref[i] <- fields[[i]][4]  ## Reference allele field
      Alt[i] <- fields[[i]][5]  ## Alternate allele field
    }
  }
  
  for(i in 1:length(Ref)) {
    ## Do a grepl match to match each possible point mutation, increment 1 to the counters
    if(grepl("A", Ref[i], fixed = TRUE) == TRUE & grepl("C", Alt[i], fixed = TRUE) == TRUE) {
      A2C <- A2C + 1
    } else 
      if(grepl("A", Ref[i], fixed = TRUE) == TRUE & grepl("G", Alt[i], fixed = TRUE) == TRUE) {
        A2G <- A2G + 1
      } else 
        if(grepl("A", Ref[i], fixed = TRUE) == TRUE & grepl("T", Alt[i], fixed = TRUE) == TRUE) {
          A2T <- A2T + 1
        } else 
          if(grepl("C", Ref[i], fixed = TRUE) == TRUE & grepl("A", Alt[i], fixed = TRUE) == TRUE) {
            C2A <- C2A + 1
          } else 
            if(grepl("C", Ref[i], fixed = TRUE) == TRUE & grepl("G", Alt[i], fixed = TRUE) == TRUE) {
              C2G <- C2G + 1
            } else 
              if(grepl("C", Ref[i], fixed = TRUE) == TRUE & grepl("T", Alt[i], fixed = TRUE) == TRUE) {
                C2T <- C2T + 1
              } else 
                if(grepl("G", Ref[i], fixed = TRUE) == TRUE & grepl("A", Alt[i], fixed = TRUE) == TRUE) {
                  G2A <- G2A + 1
                } else 
                  if(grepl("G", Ref[i], fixed = TRUE) == TRUE & grepl("C", Alt[i], fixed = TRUE) == TRUE) {
                    G2C <- G2C + 1
                  } else 
                    if(grepl("G", Ref[i], fixed = TRUE) == TRUE & grepl("T", Alt[i], fixed = TRUE) == TRUE) {
                      G2T <- G2T + 1
                    } else 
                      if(grepl("T", Ref[i], fixed = TRUE) == TRUE & grepl("A", Alt[i], fixed = TRUE) == TRUE) {
                        T2A <- T2A + 1
                      } else 
                        if(grepl("T", Ref[i], fixed = TRUE) == TRUE & grepl("C", Alt[i], fixed = TRUE) == TRUE) {
                          T2C <- T2C + 1
                        } else 
                          if(grepl("T", Ref[i], fixed = TRUE) == TRUE & grepl("G", Alt[i], fixed = TRUE) == TRUE) {
                            T2G <- T2G + 1
                          } 
                        }
  
  ## Take the sum of all point mutations
  Total_Point_Mutations <- as.numeric(sum(A2C, A2G, A2T, C2A, C2G, C2T, G2A, G2C, G2T, T2A, T2C, T2G))

  ## Calculate Percent of Each Mutation
  Mutation_Freqs <- vector(mode = "list", length = 12)
  names(Mutation_Freqs) <- c("A2C", "A2G", "A2T", "C2A", "C2G", "C2T", "G2A", "G2C", "G2T", "T2A", "T2C", "T2G")
  for(i in 1:12) {
    ## Fun trick to convert list names to the variables they represent, holding total mutation counts
    ## Then save each percent mutation to the list
    Mutation_Freqs[[i]] <- as.numeric(eval(parse(text = names(Mutation_Freqs[i])))) / Total_Point_Mutations
  }
  
  return(Mutation_Freqs)  
  
}