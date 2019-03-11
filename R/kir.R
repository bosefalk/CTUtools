#' Determine KIR Gene Copy Number
#'
#' @param kir_string String with KIR info with the first field only (three characters), example "004+003|006+010"
#' @param return_numeric If set to TRUE, returns numeric GCN instead of character, sets "x" to NA and 2|3 to 2.5 etc
#'
#' @return Default a string with "x" if "POS", "0" if NEG, NA if NA or empty string, and otherwise the GCN. If \code{return_numeric = TRUE}
#' instead returns a numeric value, where x results are NA and "2|3" are 2.5
#'
#' @examples
#'
#' str(KIR_det_GCN("004+003|006+010"))
#' str(KIR_det_GCN("004+003|006+010", return_numeric = TRUE))
#' KIR_det_GCN("POS")
#' KIR_det_GCN("POS", return_numeric = TRUE)
KIR_det_GCN = function(kir_string, return_numeric = FALSE){
  require(dplyr)

  # Call this on the KIR string if it's not empty/NA/POS/NEG:
  det_clean_GCN <- function(y) {
    split = unlist(strsplit(y,split = "[|]"))
    clean_gcn = nchar(split) - nchar(gsub("[+]","",split))+1
    clean_gcn = ifelse(min(clean_gcn) == max(clean_gcn),min(clean_gcn),paste(min(clean_gcn),max(clean_gcn),sep="|"))
    return(as.character(clean_gcn))
  }

  # Checks conditions according to top first, then working down and applying GCN determination function on string if none of
  # the other conditions are fulfilled
  gcn <- dplyr::case_when(
    grepl("POS", kir_string) ~ "x",
    grepl("NEG", kir_string) ~ "0",
    is.na(kir_string) | grepl("^$", kir_string) ~ NA_character_,
    TRUE ~ det_clean_GCN(kir_string)
  )

  if (return_numeric == TRUE) {
    # When running src/KIR_investigation.Rmd from 18-01 project error is when hitting gcn==1|2 statement no TRUE/FALSE value
    try({
      if(is.na(gcn)) {gcn <- NA_integer_}
      if(gcn == "1|2") {gcn <- 1.5}
      if(gcn == "2|3") {gcn <- 2.5}
      if(gcn == "3|4") {gcn <- 3.5}
      gcn <- suppressWarnings(as.numeric(gcn))}, silent = TRUE)
    return(gcn)
  } else {return(gcn)}
}


#' Classifies KIR 3DL1 & 3DS1 alleles into KIR3DL1-H, -L, -N
#' 
#' Uses the categorization from \code{\link{ASSIGN_KIR3DL1}}
#' 
#' @param KIR3DL1_string KIR3DL1 string with first fields only, i.e. "001/095+004". Can contain "NEG" and "POS".
#' @param KIR3DS1_string KIR3DS1 string in same format
#'
#' @return String with "KIR3DL1-H", "KIR3DL1-L", "KIR3DL1-N" or "unknown"
#' 
#' @seealso \code{\link{ASSIGN_KIR3DL1}}
#'
#' @examples
#'
#' KIR3DL1_3DS1_assignment("001/095+004", "NEG")
#' KIR3DL1_3DS1_assignment("POS", "NEG")
#' KIR3DL1_3DS1_assignment("005", "013/107")

KIR3DL1_3DS1_assignment <- function(KIR3DL1_string, KIR3DS1_string) {
  
  # This whole function was written by Henning for reading in LSL data in 1701 project. It's just inserted here pretty 
  # much as-is, with just minor tweaks to inputs and outputs (hence why the comments are in German).
  # Requires that data/ASSIGN_KIR3DL1.RData has been loaded, should have been done when loading package
  
  if (!exists("ASSIGN_KIR3DL1")) {
    load("data/ASSIGN_KIR3DL1.RData")
  }
  
  
  # Sub-functions -----------------------------------------------------------
  
  KIR_class = function(kir_glstring){
    kir_result = Extract_GLstring_KIR(kir_glstring)
    return(KIR_result_class(kir_result))
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
      result_phase = plyr::ddply(result, .(Phase), function(z){
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
  if (grepl("NEW",KIR3DS1_string) | grepl("POS",KIR3DS1_string)){
    x2 = gsub("NEW","sss",KIR3DS1_string)
    x2 = gsub("POS","sss",x2)
    kir3dl1.3ds1_GLstring = paste(KIR3DL1_string,gsub("NEW","sss",x2), sep = "+")
  }
  
  kir3dl1.3ds1_result  = KIR_class(kir3dl1.3ds1_GLstring)
  
  kir3dl1.3ds1_diff    = kir3dl1.3ds1_result[[1]]
  kir3dl1.3ds1_class   = kir3dl1.3ds1_result[[2]]
  kir3dl1.3ds1_overall = kir3dl1.3ds1_result[[3]]
  
  return(kir3dl1.3ds1_overall)
  
}
