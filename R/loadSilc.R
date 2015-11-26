loadSILC <- function(file = NULL,
                     filed = NULL,
                     filer = NULL,
                     filep = NULL,
                     fileh = NULL,
                     year = 2013,
                     country = "Austria"){
  DB135 <- NULL
  if(is.null(file)) comment1 <- "unmerged original files imported"
  if(is.null(filep)) comment1 <- "already merged original files imported"
  if(!is.null(file) & !is.null(filed)) stop("provide a merged file or unmerged files but not both")

  ##########################################
  ## read orginal file (if exists)
  if(!is.null(file)){
    nc <- nchar(file)
    if(substr(file, nc - 4, nc) == "RData" | substr(file, nc - 2, nc) == "rda"){
      xyz <- load(file)
      x <- get(xyz)
    }
    if(substr(file, nc - 2, nc) == "sav"){
      x <- haven::read_spss(file)
    }
    if(substr(file, nc - 2, nc) == "csv"){
      x <- read.csv2(file, header=TRUE, sep=",")
    }
    if(country == "France"){
      x <- subset(x, DB135 == 1)
    }
    names(x) <- tolower(names(x))
  } else {
    ##########################################
      ## read original files (if exists)
      nc <- nchar(filep)
      if(substr(filep, nc - 4, nc) == "RData" | substr(filep, nc - 2, nc) == "rda"){
          xyzp <- load(filep)
          x_p <- get(xyzp)
          xyzp <- load(fileh)
          x_h <- get(xyzp)
          xyzp <- load(filer)
          x_r <- get(xyzp)
          xyzp <- load(filed)
          x_d <- get(xyzp)
      }
      if(substr(filep, nc - 2, nc) == "sav"){
        x_p <- haven::read_spss(filep)
        x_h <- haven::read_spss(fileh)
        x_r <- haven::read_spss(filer)
        x_d <- haven::read_spss(filed)
      }
      if(substr(filep, nc - 2, nc) == "csv"){
        x_p <- read.csv2(filep, header=TRUE, sep=",")
        x_h <- read.csv2(fileh, header=TRUE, sep=",")
        x_r <- read.csv2(filer, header=TRUE, sep=",")
        x_d <- read.csv2(filed, header=TRUE, sep=",")
      }
      ## ensure that colnames are lower case
      names(x_d) <- tolower(names(x_d))
      names(x_h) <- tolower(names(x_h))
      names(x_p) <- tolower(names(x_p))
      names(x_r) <- tolower(names(x_r))

      x <- list("p"=x_p,
               "h"=x_h,
               "r"=x_r,
               "d"=x_d
     )
  }

  return(x)
}

#
# ### Test load different formats
# orig_file <- "~/workspace/sga_sdc/silcPop/new_workfile.RData"
# x <- loadSILC(orig_file)
#
# orig_file2 <- "~/workspace/sga_sdc/silcPop/new_workfile.sav"
# x <- loadSILC(orig_file)
#
# filed <- "~/workspace/sga_sdc/eurostatFormatData/2013/zielvar_d_eurostat2013.sav"
# filer <- "~/workspace/sga_sdc/eurostatFormatData/2013/zielvar_r_eurostat2013.sav"
# filep <- "~/workspace/sga_sdc/eurostatFormatData/2013/zielvar_p_eurostat2013.sav"
# fileh <- "~/workspace/sga_sdc/eurostatFormatData/2013/zielvar_h_eurostat2013.sav"
# suf4 <- loadSILC(filed = filed, filer = filer, filep = filep, fileh = fileh)

mergeSILC <- function(d, r, h, p){
  ## x ... list with four elements
  pers <- merge(x=r, y=p, by.x="rb030", by.y="pb030", all.x=TRUE)
  pers$rb030 <- pers$pb030
  hh <- merge(x=d, y=h, by.x="db030", by.y="hb030", all.x=TRUE)
  hh$hb030 <- hh$db030
  if(!("rx030" %in% colnames(pers))) pers$rx030 <- pers$rb040
  if(!("rx030" %in% colnames(hh))) hh$rb040 <- hh$pb030
  if(!("hb030" %in% colnames(hh))) hh$hb030 <- pers$db030
  m <- merge(x=pers, y=hh, by.x="rx030", by.y="hb030")
  return(m)
}

# suf <- mergeSILC(suf4[["d"]], suf4[["r"]], suf4[["h"]], suf4[["p"]])


checkCol <- function(x, y){
  int <- intersect(colnames(x), colnames(y))
  w <- which(!(colnames(y) %in% colnames(x)))
  cn <- colnames(y)[w]
  cat("\n these variables are in y but not in x:\n", cn, "\n")
  return(cn)
}

# checkCol(x, suf)


chooseSILCvars <- function(x, vars = c("db030", "db040", "rb030", "rb080", "rb090", "pl031", "pb220a", "py010g",
                                       "py021g", "py050g",  "py080g", "py090g", "py100g", "py110g", "py120g", "py130g",
                                       "py140g", "hy040g", "hy050g", "hy060g",  "hy070g", "hy080g" , "hy090g" ,"hy100g" ,
                                       "hy110g", "hy120g", "hy130g", "hy140g", "db090", "rb050", "pb190", "pe040",
                                       "pl051","pl111" , "rb010"), country = NULL){
  x <- x[, vars]
  if(!is.numeric(x$db090)) x$db090 <- as.numeric(as.character(x$db090))
  if(!is.numeric(x$rb050)) x$rb050 <- as.numeric(as.character(x$rb050))
  if(!is.factor(x$pb190)) x$pb190 <- factor(x$pb190)
  if(!is.factor(x$pb190)) x$pb190 <- factor(x$pb190)
  if(!is.factor(x$pe040)) x$pe040 <- factor(x$pe040)
  if(!is.factor(x$pl051)) x$pl051 <- factor(x$pl051)
  if(!is.factor(x$pl111)) x$pl111 <- factor(x$pl111)
  if(!is.factor(x$rb010)) x$rb010 <- factor(x$rb010)
  if(!is.factor(x$pl031)) x$pl031 <- factor(x$pl031)

  ## category 1 is too small:
  tab <- table(x$pe040, useNA = "always")[2]
  if(tab < 10){
    x$pe040 <- plyr::revalue(x$pe040, c("0"="0-1", "1"="0-1" ))
    message(cat("Note: number of categories in pe040 was too small ( count =", tab, ").\n Categories 0 and 1 in pe040 have been combined to category 0-1\n"))
  }

  if(country == "Austria"){
    x$db040 <- factor(x$db040,
                           labels= c("Burgenland", "Lower Austria", "Vienna",
                                     "Carinthia", "Styria", "Upper Austria",
                                     "Salzburg", "Tyrol", "Vorarlberg"))
  }

  x <- x[order(x$db030),]
  x$db030 <- restructureHHid(x)
  return(x)
}


modifySILC <- function(x, country = "Austria"){
  x$age <- getAge(data = x)
  x$rb090 <- getGender(data = x)
  x$hsize <- getHsize(data = x)
  x$hsize <- as.factor(x$hsize)
  x$pl031 <- getEcoStat(data = x, levels = levels(x$pl031))
  ## modify pb220a
  if(country == "France"){
    owncountry <- "FR"
    EU <- c("AT","BE","BG","CY","CZ", "DE", "DK","EE","EL","ES","FI","GR","HU","IE",
            "IT","LT","LU","LV","MT","NL","PL","PT","RO","SI","SE","SK","UK")
  }
  if(country == "Austria"){
    owncountry <- "AT"
    EU <- c("BE","BG","CY","CZ", "DE", "DK","EE","EL","ES","FI","FR","GR","HU","IE",
            "IT","LT","LU","LV","MT","NL","PL","PT","RO","SI","SE","SK","UK")
  }
  other <- c("CAN","CH","CSA","HR","IS","ME","MK","NAF","NME","NO",
             "OAF","OAS","OCE","OEU","OT","OTH","TR","USA","WAF")
  x$pb220a <- factor(x$pb220a)
  x$pb220a <- getCitizenship(data=x, owncountry=owncountry,
                                  EU=EU, other=other)

  ## Recoding of Occupation to 1-digit version
  x$pl051 <- as.numeric(as.character(x$pl051))
  x$pl051 <- as.factor(trunc(x$pl051/10))
  x$pl051 <- plyr::revalue(x$pl051, c("0"="0-1", "1"="0-1" ))

  if(country %in% c("France", "Austria")){
    ## Recoding of NACE code - this may be country-specific.
    levels(x$pl111)[which(levels(x$pl111) %in% c("1","2","3"))] <- "a"
    levels(x$pl111)[which(levels(x$pl111) %in% c("5","6","7","8","9",
                                                           "10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25",
                                                           "26","27","28","29","30","31","32","33","34","35","36","37","38","39"))] <- "b-e"
    levels(x$pl111)[which(levels(x$pl111) %in% c("41","42","43"))] <- "f"
    levels(x$pl111)[which(levels(x$pl111) %in% c("45","46","47"))] <- "g"
    levels(x$pl111)[which(levels(x$pl111) %in% c("49","50","51","52",
                                                           "53"))] <- "h"
    levels(x$pl111)[which(levels(x$pl111) %in% c("55","56"))] <- "i"
    levels(x$pl111)[which(levels(x$pl111) %in% c("58","59","60","61",
                                                           "62","63"))] <- "j"
    levels(x$pl111)[which(levels(x$pl111) %in% c("64","65","66"))] <- "k"
    levels(x$pl111)[which(levels(x$pl111) %in% c("68","69","70",
                                                           "71","72","73","74","75","76","77","78","79","80","81","82"))] <- "l-n"
    levels(x$pl111)[which(levels(x$pl111)=="84")] <- "o"
    levels(x$pl111)[which(levels(x$pl111)=="85")] <- "p"
    levels(x$pl111)[which(levels(x$pl111) %in% c("86","87","88"))] <- "q"
    levels(x$pl111)[which(levels(x$pl111) %in% c("90","91","92","93","94",
                                                           "95","96","97","98","99"))] <- "r-u"
  }
  if(country == "Germany"){
    # Code for German recoding:
    levels(x$pl111)[which(levels(x$pl111)=="1")] <- "a"
    levels(x$pl111)[which(levels(x$pl111) %in% c("2","3","4","5", "6"))] <- "b-e"
    levels(x$pl111)[which(levels(x$pl111)=="7")] <- "f"
    levels(x$pl111)[which(levels(x$pl111)=="8")] <- "g"
    levels(x$pl111)[which(levels(x$pl111)=="9")] <- "h"
    levels(x$pl111)[which(levels(x$pl111)=="10")] <- "i"
    levels(x$pl111)[which(levels(x$pl111)=="11")] <- "j"
    levels(x$pl111)[which(levels(x$pl111)=="12")] <- "k"
    levels(x$pl111)[which(levels(x$pl111) %in% c("13","14","15"))] <- "l-n"
    levels(x$pl111)[which(levels(x$pl111)=="16")] <- "o"
    levels(x$pl111)[which(levels(x$pl111)=="17")] <- "p"
    levels(x$pl111)[which(levels(x$pl111)=="18")] <- "q"
    levels(x$pl111)[which(levels(x$pl111) %in% c("19","20","21", "22","23"))] <- "r-u"
}

  x$pgrossIncome <- rowSums(x[, c("py010g","py021g","py050g","py080g",
                                            "py090g","py100g", "py110g","py120g", "py130g","py140g")], na.rm = TRUE)
  x$hgrossIncome <- rowSums(x[, c("hy040g","hy050g","hy060g",
                                            "hy070g", "hy080g", "hy090g", "hy110g")], na.rm = TRUE)
  x$hgrossminus <- rowSums(x[, c("hy120g","hy130g","hy140g")], na.rm = TRUE)

  breaks <- getBreaks(x$pgrossIncome, x$rb050,
                      upper = Inf, equidist = FALSE, zeros = TRUE)
  x$pgrossIncomeCat <- getCat(x$pgrossIncome, breaks)
  breakshh <- getBreaks(x$hgrossIncome, x$rb050,
                        upper = Inf, equidist = FALSE, zeros = TRUE)
  x$hgrossIncomeCat <- getCat(x$hgrossIncome, breakshh)
  breakshhm <- getBreaks(x$hgrossminus, x$rb050,
                         upper = Inf, equidist = FALSE, zeros = TRUE)
  x$hgrossminusCat <- getCat(x$hgrossminus, breakshhm)

  x$country <- factor(rep(country, nrow(x)))

  if(country == "France"){
    x$db040[which(x$db030 %in% c(9382,11461))] <- "FR82"
    x$db040 <- droplevels(x$db040)
  }

  return(x)
}

# x2 <- chooseSILCvars(x, country = "Austria")
# dim(x2)
# x2 <- modifySILC(x2)
# dim(x2)




