#' Determine KIR Gene Copy Number
#'
#' @param kir_string String with KIR info with the first field only (three characters), example "004+003|006+010", "NEG", "POS"
#' @param return_numeric If set to TRUE, returns numeric GCN instead of character, sets "POS" to NA and 2|3 to 2.5 etc
#'
#' @return String vector with the gene copy number, with "x" if "POS", "0" if NEG, NA if NA or empty string. If \code{return_numeric = TRUE}
#' instead returns a numeric vector, where POS results are NA and "2|3" are 2.5
#'
#' @examples
#'
#' str(KIR_det_GCN("004+003|006+010"))
#' str(KIR_det_GCN("004+003|006+010", return_numeric = TRUE))
#' KIR_det_GCN("POS")
#' KIR_det_GCN("POS", return_numeric = TRUE)
#' 
#' dat <- data.frame(kircol = c("004+003", "001", "008+008+008", "POS", "NEG", "", NA, "001|001+002"))
#' KIR_det_GCN(dat$kircol)
#' KIR_det_GCN(dat$kircol, return_numeric = TRUE)
KIR_det_GCN = function(kir_string, return_numeric = FALSE){
  require(dplyr)

  kir_string <- as.character(kir_string)
  
  # Call this on each KIR string
  det_clean_GCN <- function(y) {
    if (is.na(y) | grepl("^$", y)) return(NA_character_)
    if (grepl("POS", y)) return("x")
    if (grepl("NEG", y)) return("0")
    
    split = unlist(strsplit(y,split = "[|]"))
    clean_gcn = nchar(split) - nchar(gsub("[+]","",split))+1
    clean_gcn = ifelse(min(clean_gcn) == max(clean_gcn),min(clean_gcn),paste(min(clean_gcn),max(clean_gcn),sep="|"))
    return(as.character(clean_gcn))
  }

  gcn <- unlist(lapply(X = kir_string, FUN = det_clean_GCN))
  

  if (return_numeric == TRUE) {
    # When running src/KIR_investigation.Rmd from 18-01 project error is when hitting gcn==1|2 statement no TRUE/FALSE value
    gcn_num <- function(y) {
      if(is.na(y)) return(NA_real_)
      if(y == "1|2") return(1.5)
      if(y == "2|3") return(2.5)
      if(y == "3|4") return(3.5)
      return(suppressWarnings(as.numeric(y)))
    }
    gcn <- unlist(lapply(X = gcn, FUN = gcn_num))
    return(gcn)
  } else {return(gcn)}
}


#' Classifies KIR 3DL1 & 3DS1 alleles into KIR3DL1-H, -L, -N
#' 
#' Uses the categorization from \code{\link{ASSIGN_KIR3DL1}}
#' 
#' @param KIR3DL1_string KIR3DL1 string vector with first fields only, i.e. "001/095+004". Can contain "NEG" and "POS".
#' @param KIR3DS1_string KIR3DS1 string vector in same format
#'
#' @return String vector with "KIR3DL1-H", "KIR3DL1-L", "KIR3DL1-N" or "unknown"
#' 
#' @seealso \code{\link{ASSIGN_KIR3DL1}}
#'
#' @examples
#'
#' dat <- data.frame(kir_3DL1 = c("001/095+004", "POS", "005"), kir_3DS1 = c("NEG", "NEG", "013/107"), stringsAsFactors = FALSE)
#' KIR3DL1_3DS1_assignment(dat$kir_3DL1, dat$kir_3DS1)
KIR3DL1_3DS1_assignment <- function(KIR3DL1_string, KIR3DS1_string) {
  
  # If passed as factor variables pass to character
  KIR3DL1_string <- as.character(KIR3DL1_string)
  KIR3DS1_string <- as.character(KIR3DS1_string)
  
  # This whole function was written by Henning for reading in LSL data in 1701 project. It's just inserted here pretty 
  # much as-is, with just minor tweaks to inputs and outputs (hence why the comments are in German).
  # Requires that data/ASSIGN_KIR3DL1.RData has been loaded, should have been done when loading package
  
  if (!exists("ASSIGN_KIR3DL1")) {
    load("data/ASSIGN_KIR3DL1.RData")
  }
  
  
  
  # Sub-functions -----------------------------------------------------------
  
  KIR_class = function(kir_glstring){
    kir_result = lapply(kir_glstring, Extract_GLstring_KIR)
    #kir_result <- Extract_GLstring_KIR(kir_glstring)
    #KIR_result_class(kir_result)
    out = lapply(kir_result, KIR_result_class)
    return(out)
  }
  
  
  ## Aufl?sen des GL-Strings
  Extract_GLstring_KIR = function(kir_glstring){
    # 1. nach | Phasings auftrennen,
    kir_result = {}
    kir_glstring = gsub("[*]","",kir_glstring) # "*"-Zeichen l?schen, wenn vorhanden
    x1 = unlist(strsplit(kir_glstring, split = "|", fixed = T))
    for (i in 1:length(x1)){
      # 2. nach + beiden Allelen auftrennen
      x2 = unlist(strsplit(x1[i], split = "+", fixed = T))
      for (j in 1:length(x2)){
        # 3. nach / verschiedenen Allele auftrennen
        x3 = unlist(strsplit(x2[j], split = "/", fixed = T))
        #paste(unique(substr(x3,start = 1, stop = 3)), collapse = "/")
        x3 = unique(substr(x3,start = 1, stop = 3))
        for (k in 1:length(x3)) kir_result = rbind(kir_result,c(i,j,x3[k]))
      }
    }
    kir_result = as.data.frame(kir_result)
    names(kir_result) = c("Phase","Allel","Digit")
    kir_result = merge(kir_result, ASSIGN_KIR3DL1, all.x = T)
    kir_result$Class[is.na(kir_result$Class)] = 0
    # kir_result$Class = sapply(kir_result$Digit,function(x){
    #   y = ASSIGN_KIR3DL1$Class[ASSIGN_KIR3DL1$Allele==x]
    #   if (length(y) == 0) y = 0
    #   # if (y == 0 & x !="NEW") print(x) # Ausgabe, wenn ein neues Allel gefunden wird
    #   return(y)
    # })
    return(kir_result)
  }
  
  
  
  ## Zuordnungstabelle f?r die KIR3DL1/KIR3DS1-Klassifizierung
  assignment_KIR3DL1 = data.frame(subtype = c("l","h","s","n"), 
                                  num = 1:4, 
                                  overall = c("KIR3DL1-L","KIR3DL1-H",rep("KIR3DL1-N",2)))
  
  KIR_result_class = function(kir_result){
    xx = subset(kir_result, select = - Digit)
    xx = xx[!duplicated(xx),] #dadurch werden mehrdeutige Allele, die aber diesselbe Subtyp-Klasse l,h,s oder n haben, zusammengefasst  
    result = plyr::ddply(xx,c("Phase","Allel"),function(z){
      c(subtype = paste(z$Class,collapse = "/"))
    })
    output = plyr::ddply(result,c("Phase"),function(z){
      paste(z$subtype,collapse = "+")
    })
    output = paste(output$V1,collapse = "|")  
    
    
    if (any(result$subtype == 0) | any(grepl("/",result$subtype))){
      if (any(result$subtype == 0)){
        class = "unclassified allel"
        overall = "unknown"
      }
      if (any(grepl("/",result$subtype))){
        class = "ambiguites"
        overall = "unknown"
      }
    } else {
      result$Phase = as.numeric(result$Phase)
      ## !!! Annahme: Bei 1 Allel wird das zweite zu "n" gesetzt
      for (i in 1:max(result$Phase)){
        if (length(result$Phase[result$Phase == i]) == 1){
          result = rbind(result,c(i,2,"n"))
        }
      }
      result = plyr::join(result,assignment_KIR3DL1,by = "subtype") 
      ## !!! Annahme: Bei mind. 3 Allelen werden die beiden st?rksten Inhibierungseffekte gesucht
      result_phase = plyr::ddply(result, "Phase", function(z){
        z = z[order(z$num)[1:2],]
        c(class_i = as.character(paste(z$subtype,collapse = "/")),
          overall_i = as.character(z$overall[1]))
      })
      class = ifelse(length(unique(result_phase$class_i))==1,result_phase$class_i[1],"phasing")
      overall = ifelse(length(unique(result_phase$overall_i))==1,result_phase$overall_i[1],"unknown") 
    }
    return(list(as.character(output),as.character(class),as.character(overall)))
  }
  
  # main call ---------------------------------------------------------------
  
  
  kir3dl1.3ds1 = ifelse(KIR3DS1_string == "NEG", KIR3DL1_string, paste(KIR3DL1_string,KIR3DS1_string, sep = "+"))
  kir3dl1.3ds1_GLstring = kir3dl1.3ds1
  
  
  
  ## Wenn kir3ds1 positiv ist oder ein neues Allel typisiert wurde, wird per se ein "s" klassifiziert
  kir3dl1.3ds1_GLstring <- unlist(lapply(kir3dl1.3ds1_GLstring, FUN = function (x) {
    if (grepl("NEW",x) | grepl("POS",x)){
    x2 = gsub("NEW","sss",x)
    x2 = gsub("POS","sss",x2)
    out = paste(x,gsub("NEW","sss",x2), sep = "+")
    return(out)
    } else {
    return(x)
  }
    }))
  
  kir3dl1.3ds1_result  = KIR_class(kir3dl1.3ds1_GLstring)
  
  #kir3dl1.3ds1_diff    = kir3dl1.3ds1_result[[1]]
  #kir3dl1.3ds1_class   = kir3dl1.3ds1_result[[2]]
  #kir3dl1.3ds1_overall = kir3dl1.3ds1_result[[3]]
  
  overall <- unlist(lapply(kir3dl1.3ds1_result, function(x) return(x[[3]])))
  
  return(overall)
  
}

#' Classifies KIR3DL1 and HLA-B into strong/weak/no inhibiting
#' 
#' 
#' @param KIRD3L1_assignment String vector, KIR3DL1 assignment in fromat "KIR3DL1-H", from \code{\link{KIR3DL1_3DS1_assignment}}
#' @param HLA_B_overall String vector with HLA-B Bw classification in format "Bw4 - 80T", from \code{\link{HLA_B_classification}}
#' @param levels Can be set to 2, 3 or 4, default is 2
#'
#' @return If levels = 2, string with "Strong inhibiting", "Weak inhibiting/noninhibiting" or "unknown". If levels = 3, 
#' split into "Weak inhibiting" and "non-inhibiting". If levels = 4, split "noninhibiting" into "Missing ligand" and "Educated, uninhibiting"
#' 
#'
#' @examples
#' dat <- data.frame(assignment_KIR3DL1 = c("KIR3DL1-L", "KIR3DL1-N", "unknown", "unknown"), HLA_B_group = c("Bw4 - 80T", "Bw4 - 80I", NA, "Bw4 - 80I"), stringsAsFactors = FALSE)
#' KIR3DL1_HLA_B_inhibiting(dat$assignment_KIR3DL1, dat$HLA_B_group)
#' KIR3DL1_HLA_B_inhibiting(dat$assignment_KIR3DL1, dat$HLA_B_group, levels = 3)
#' KIR3DL1_HLA_B_inhibiting(dat$assignment_KIR3DL1, dat$HLA_B_group, levels = 4)
KIR3DL1_HLA_B_inhibiting <- function(KIR3DL1_assignment, HLA_B_overall, levels = 2) {
  # If passed as factor columns
  KIR3DL1_assignment <- as.character(KIR3DL1_assignment)
  HLA_B_overall <- as.character(HLA_B_overall)
  
  # dplyr::case_when used instead of lots of ifelse statements. If a left hand side evaluates to TRUE, the right hand side
  # value is returned, otherwise it continues down the list. 
  
  # First constructs the most detailed case, when levels = 4. Then, if levels is set lower it merges groups together.
    out <- case_when(
      is.na(KIR3DL1_assignment) | is.na(HLA_B_overall) ~ NA_character_,
      KIR3DL1_assignment == "unknown" | HLA_B_overall == "unknown" ~ "unknown",
      KIR3DL1_assignment == "KIR3DL1-L" & HLA_B_overall == "Bw4 - 80T" ~ "Strong inhibiting",
      KIR3DL1_assignment == "KIR3DL1-H" & HLA_B_overall == "Bw4 - 80I" ~ "Strong inhibiting",
      KIR3DL1_assignment == "KIR3DL1-L" & HLA_B_overall == "Bw4 - 80I" ~ "Weak inhibiting",
      KIR3DL1_assignment == "KIR3DL1-H" & HLA_B_overall == "Bw4 - 80T" ~ "Weak inhibiting",
      KIR3DL1_assignment == "KIR3DL1-L" & HLA_B_overall == "Bw6" ~ "Missing ligand",
      KIR3DL1_assignment == "KIR3DL1-H" & HLA_B_overall == "Bw6" ~ "Missing ligand",
      KIR3DL1_assignment == "KIR3DL1-N" & HLA_B_overall == "Bw6" ~ "Missing ligand",
      KIR3DL1_assignment == "KIR3DL1-N" & HLA_B_overall == "Bw4 - 80I" ~ "Educated, Uninhibited",
      KIR3DL1_assignment == "KIR3DL1-N" & HLA_B_overall == "Bw4 - 80T" ~ "Educated, Uninhibited"
    )

    if (levels == 4) {return(out)}
    
    # group Missing Ligand and Educated, Uninhibited into "noninhibiting"
      out <- ifelse((out == "Missing ligand" | out == "Educated, Uninhibited"),
             "noninhibiting", out)
    
    if (levels == 3) {return(out)}
    
    # For levels = 2, also group weak and noninhibiting together
      out <- ifelse((out == "noninhibiting" | out == "Weak inhibiting"),
                    "Weak inhibiting/noninhibiting", out)
    return(out) # Returns here if levels = 2
}



#' KIR present/absent
#'
#' Determine whether the KIR is present, given a KIR string. Assumes a string with only "NEG" is absent, any other string not including
#' "NEG" is present, if missing or some combination together with "NEG" return NA - see examples.
#'
#' @param string KIR string vector
#'
#' @return logical vector, TRUE if the KIR is present, FALSE if it's not present, NA if passed an empty/NA string
#'
#' @examples
#' dat <- data.frame(kirstring = c("003/034+003/034", "NEG", "", NA, "003|NEG"), stringsAsFactors = FALSE)
#' KIR_present(dat$kirstring)
KIR_present <- function(string) {
  # if passed factor column
  string <- as.character(string)
  # Assumes anything other than an empty string or a "NEG" KIR string means the gene is present
  # Need to make sure inputs are cleaned
  out <- case_when(
    is.na(string) ~ NA,
    grepl("^$", string) ~ NA,
    string == "NA" ~ NA,
    grepl("^NEG$", string) ~ FALSE, # Only NEG
    grepl("NEG", string) ~ NA, # NEG with something else
    TRUE ~ TRUE) # Otherwise present
  return(out)
}




#' KIR Centromeric Haplotype
#' 
#' Assigns the Centrometic haplotype (A/A, A/B or B/B) given a set of KIRs
#'
#' @param df data.frame with KIR alleles string as columns
#' @param KIR_2DL3col string, name of 2DL3 column
#' @param KIR_2DS2col string, name of 2DS2 column
#' @param KIR_2DL2col string, name of 2DL2 column
#' @param present_cols Default is FALSE which assumes the KIR columns contains the allele string ("001", "NEG" etc), and the presence / absence is calculated inside the function using \code{\link{KIR_present}}. If set to TRUE, the KIR columns already contain TRUE/FALSE for presence/absence for all KIRs. 
#'
#' @return string vector with haplotype ("A/A", "A/B", "B/B"), NA if any input KIR is missing, and "unknown" otherwise.
#'
#' @seealso \code{\link{KIR_haplotype_Tel}}
#' @examples
#' dat <- data.frame(kir_2DL3 = c("001", "NEG", NA, "001|002+003"), 
#'    kir_2DS2 = c("001", "001", "001", "001"), 
#'    kir_2DL2 = c("NEG", "001", "NEG", "NEG"), 
#'    stringsAsFactors = FALSE)
#' KIR_haplotype_Cen(dat, "kir_2DL3", "kir_2DS2", "kir_2DL2")
#' 
#' dat_pres <- data.frame(kir_2DL3 = c(TRUE, FALSE, NA, TRUE), 
#'    kir_2DS2 = c(TRUE, TRUE, TRUE, TRUE), 
#'    kir_2DL2 = c(FALSE, TRUE, FALSE, FALSE), 
#'    stringsAsFactors = FALSE)
#' KIR_haplotype_Cen(dat_pres, "kir_2DL3", "kir_2DS2", "kir_2DL2", present_cols = TRUE)
KIR_haplotype_Cen <- function(df, KIR_2DL3col, KIR_2DS2col, KIR_2DL2col, present_cols = FALSE) {
  
  if (present_cols == FALSE) {
    if (any(class(df[[KIR_2DL3col]]) == "logical", 
            class(df[[KIR_2DS2col]]) == "logical",
            class(df[[KIR_2DL2col]]) == "logical")) {
      stop("KIR columns are not 'logical' type but present_cols was set to TRUE")
    }
    kir2DL3_present <- CTUtools::KIR_present(as.character(df[[KIR_2DL3col]]))
    kir2DS2_present <- CTUtools::KIR_present(as.character(df[[KIR_2DS2col]]))
    kir2DL2_present <- CTUtools::KIR_present(as.character(df[[KIR_2DL2col]]))
  }
  if (present_cols == TRUE) {
    if (any(class(df[[KIR_2DL3col]]) != "logical", 
            class(df[[KIR_2DS2col]]) != "logical",
            class(df[[KIR_2DL2col]]) != "logical")) {
      stop("KIR columns are not 'logical' type but present_cols was set to TRUE")
    }
    kir2DL3_present <- df[[KIR_2DL3col]]
    kir2DS2_present <- df[[KIR_2DS2col]]
    kir2DL2_present <- df[[KIR_2DL2col]]
  }
  
  out <- case_when(
    is.na(kir2DL3_present) | is.na(kir2DS2_present) | is.na(kir2DL2_present) ~ NA_character_,
    kir2DL3_present == "" | kir2DS2_present == "" | kir2DL2_present == "" ~ NA_character_,
    (kir2DL3_present == TRUE & kir2DS2_present == FALSE & kir2DL2_present == FALSE) ~ "A/A",
    (kir2DL3_present == TRUE & (kir2DS2_present == TRUE | kir2DL2_present == TRUE)) ~ "A/B",
    (kir2DL3_present == FALSE & (kir2DS2_present == TRUE | kir2DL2_present == TRUE)) ~ "B/B",
    TRUE ~ "unknown"
  )
  
  return(out)
  
}


#' KIR Telomeric Haplotype
#' 
#' Assigns the Telomeric haplotype (A/A, A/B or B/B) given a set of KIRs
#'
#' @param df data.frame with KIR alleles string as columns
#' @param KIR_3DL1col string, name of 3DL1 column
#' @param KIR_2DS4col string, name of 2DS4 column
#' @param KIR_3DS1col string, name of 3DS1 column
#' @param KIR_2DS1col string, name of 2DS1 column
#' @param present_cols Default is FALSE which assumes the KIR columns contains the allele string ("001", "NEG" etc), and the presence / absence is calculated inside the function using \code{\link{KIR_present}}. If set to TRUE, the KIR columns already contain TRUE/FALSE for presence/absence for all KIRs. 
#'
#' @return string vector with haplotype ("A/A", "A/B", "B/B"), NA if any input KIR is missing, and "unknown" otherwise.
#'
#' @seealso \code{\link{KIR_haplotype_Cen}}
#' @examples
#' dat <- data.frame(kir_3DL1 = c("001", "NEG", NA, "001|002+003"), 
#'    kir_2DS4 = c("001", "001", "001", "001"), 
#'    kir_3DS1 = c("NEG", "001", "NEG", "NEG"),
#'    kir_2DS1 = c("001", "001", "NEG", "NEG"),
#'    stringsAsFactors = FALSE)
#' KIR_haplotype_Tel(dat, "kir_3DL1", "kir_2DS4", "kir_3DS1", "kir_2DS1")
#'    
#' dat_pres <- data.frame(kir_3DL1 = c(TRUE, FALSE, NA, TRUE), 
#'     kir_2DS4 = c(TRUE, TRUE, TRUE, TRUE), 
#'     kir_3DS1 = c(FALSE, TRUE, FALSE, FALSE),
#'     kir_2DS1 = c(TRUE, TRUE, FALSE, FALSE),
#'     stringsAsFactors = FALSE)
#' KIR_haplotype_Tel(dat_pres, "kir_3DL1", "kir_2DS4", "kir_3DS1", "kir_2DS1", present_cols = TRUE)
KIR_haplotype_Tel <- function(df, KIR_3DL1col, KIR_2DS4col, KIR_3DS1col, KIR_2DS1col, present_cols = FALSE) {
  
  if (present_cols == FALSE) {
    if (any(class(df[[KIR_3DL1col]]) == "logical", 
            class(df[[KIR_2DS4col]]) == "logical",
            class(df[[KIR_3DS1col]]) == "logical",
            class(df[[KIR_2DS1col]]) == "logical")) {
      stop("KIR column are 'logical' type but present_cols was set to FALSE")
    }
    kir3DL1_present <- CTUtools::KIR_present(as.character(df[[KIR_3DL1col]]))
    kir2DS4_present <- CTUtools::KIR_present(as.character(df[[KIR_2DS4col]]))
    kir3DS1_present <- CTUtools::KIR_present(as.character(df[[KIR_3DS1col]]))
    kir2DS1_present <- CTUtools::KIR_present(as.character(df[[KIR_2DS1col]]))
  }
  if (present_cols == TRUE) {
    if (any(class(df[[KIR_3DL1col]]) != "logical", 
            class(df[[KIR_2DS4col]]) != "logical",
            class(df[[KIR_3DS1col]]) != "logical",
            class(df[[KIR_2DS1col]]) != "logical")) {
      stop("KIR columns are not 'logical' type but present_cols was set to TRUE")
    }
    kir3DL1_present <- df[[KIR_3DL1col]]
    kir2DS4_present <- df[[KIR_2DS4col]]
    kir3DS1_present <- df[[KIR_3DS1col]]
    kir2DS1_present <- df[[KIR_2DS1col]]
  }
  
  
  
  out <- case_when(
    is.na(kir3DL1_present) | is.na(kir2DS4_present) | is.na(kir3DS1_present) | is.na(kir2DS1_present) ~ NA_character_,
    kir3DL1_present == "" | kir2DS4_present == "" | kir3DS1_present == "" | kir2DS1_present == "" ~ NA_character_,
    (kir3DL1_present == TRUE & kir2DS4_present == TRUE) & (kir3DS1_present == FALSE & kir2DS1_present == FALSE) ~ "A/A",
    (kir3DL1_present == TRUE & kir2DS4_present == TRUE) & (kir3DS1_present == TRUE | kir2DS1_present == TRUE) ~ "A/B",
    (kir3DL1_present == FALSE | kir2DS4_present == FALSE) ~ "B/B",
    TRUE ~ "unknown"
  )
  
  return(out)
  
}

