#' Convert hap objects to other formats
#'
#' Takes in hap files and converts them to the formats needed for other analysis programs
#'
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#'
#' @param hap the hap object to convert
#' @param format the format you wish to convert the hap object to
#' @param data.col the column corresponding to the first individual
#' @param write.file should the converted file be written
#' @param filename name to use for the output file
#' @param parents columns for the two parents being used to convert to RQTL and AB format
#' @param encoding numbers to reencode alleles (minor, het, major)
#' @param genotypes symbols used to denote genotype calls in the hap object
#' @param missing symbols used to denote missing data in the hap object
#' @param jm.pop joinmap population type
#'
#' @keywords
#'
#' @examples
#'
#' @export

hap.convert <- function(hap, format, data.col = 13, write.file=FALSE, filename, parents=NULL, encoding=c(-1,0,1), genotypes=c("A","C","G","T","H"), missing=c("N",NA), jm.pop) {

  # TODO add: PLINK (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed)
  # TODO add: CERVUS: three files (list possible parents, list of all progeny, genotype file)

  if(missing(format)) {
    stop("No export format specified.")
  }
  
  if(write.file==TRUE & missing(filename)) {
    stop("Must specify base filename if writing the output file.")
  }

  MEGA.F = function(...) {
    # TODO http://www.megasoftware.net/mega4/mega4.pdf
    
    if(write.file) {
      write.converted(hap.calls, file.name, "mega", ".txt", sep="\t")
    }
  }

  RQTL.F = function(...) {
    # TODO http://www.inside-r.org/packages/cran/qtl/docs/read.cross

    # line id in first column, marker names across header, chromosome on second line, A B H genotypes, specific columns empty
    
    if(write.file) {
      write.converted(hap.calls, file.name, "rqtl", ".txt", sep="\t")
    }
  }

  STRUCTURE.F = function(...) {
    strfile = hap
    
    # remove SNPs where 2 or more SNPs in same tag / keeps the first SNP in that tag
    if(!"rsOrig"%in%colnames(strfile)) {
      strfile = strfile[!duplicated(strfile$rsOrig),]
    }
    
    strfile = as.matrix(strfile[,data.col:ncol(strfile)])
    
    # Transpose the strfile for STRUCTURE input
    strfile <- t(strfile)
    
    # Replace genotypes with numeric values
    for(i in 1:length(genotypes)) {
      strfile[strfile==genotypes[i]] = i
    }
    
    for(i in 1:length(missing)) {
      strfile[strfile==missing[i]] = -9
      
      if(is.na(missing[i])) {
        strfile[is.na(strfile)] = -9
      }
    }
    
    as.data.frame(strfile)
    
    # writing out the strfile for STRUCTURE program in .txt format
    if(write.file) {
      write.converted(strfile, file.name, "structure", ".txt", sep="\t")
    }
  }

  FSTRUCTURE.F = function(...) {
    # TODO https://github.com/rajanil/fastStructure

    # rows in the data file correspond to samples, with two rows per sample and columns correspond to SNPs
    # The first 6 columns of the file will be ignored; typically include IDs, metadata, etc.
    # only handles bi-allelic loci. The two alleles at each locus can be encoded as desired
    # missing data should be encoded as -9.
    
    if(write.file) {
      write.converted(hap.calls, file.name, "fstructure", ".txt", sep="\t")
    }
  }

  AB.F = function(...) {
    check.parents(parents,hap)

    # If user specifies column number instead of names
    if(is.numeric(parents)) {
      p1col = parents[1]
      p2col = parents[2]

      print(paste("Using ",colnames(hap)[parents[1]]," and ",colnames(hap)[parents[2]]," as parents.",sep=""))
    } else {
      p1col = which(grepl(parents[1],colnames(hap),ignore.case=T))
      p2col = which(grepl(parents[2],colnames(hap),ignore.case=T))
    }

    ## Remove markers that are hets, identical, or missing in both parents
    print("Removing markers that are het, identical, or missing in both parents...")

    hap.clean = hap[,data.col:ncol(hap)]
    hap.clean = hap[hap[,p1col]!="H" & hap[,p2col]!="H",]
    hap.clean = hap.clean[hap.clean[,p1col]!="N" & hap.clean[,p2col]!="N",]
    hap.clean = hap.clean[(hap.clean[,p1col]!=hap.clean[,p2col]),]

    print(paste("Removed ",nrow(hap)-nrow(hap.clean)," markers.",sep=""))

    ## Identify two alleles for each marker
    if(!"alleles"%in%colnames(hap.clean)) {
      alleles = derive.alleles(hap.clean)
      allele1 = alleles$allele1
      allele2 = alleles$allele2
    } else {
      allele1 = substring(hap$alleles,1,1)
      allele2 = substring(hap$alleles,3,3)
    }

    ## Identify which allele belongs to which parent and impute missing from the other parent
    hap.ab = hap.clean

    hap.ab[,p1col] = as.character(hap.ab[,p1col])
    hap.ab[,p2col] = as.character(hap.ab[,p2col])

    print("Imputing parental genotypes when only one present...")

    for(i in 1:length(alleles)){
      if (hap.ab[,p1col][i] == "N") {
        if (hap.ab[,p2col][i] == allele1[i]) {
          hap.ab[,p1col][i] = allele2[i]
        } else if (hap.ab[,p2col][i] == allele2[i]) {
          hap.ab[,p2col][i] = allele1[i]
        }
      } else if (hap.ab[,p2col][i] == "N") {
        if (hap.ab[,p1col][i] == allele1[i]) {
          hap.ab[,p2col][i] = allele2[i]
        } else if (hap.ab[,p1col][i] == allele2[i]) {
          hap.ab[,p2col][i] = allele1[i]
        }
      }
    }

    ## Convert to A/B
    hap.calls = hap.ab
    hap.calls <- data.frame(lapply(hap.calls, as.character), stringsAsFactors=FALSE, check.names=FALSE)
    hap.calls[hap.calls==hap.ab[,p2col]] = "B"
    hap.calls[hap.calls==hap.ab[,p1col]] = "A"

    hap.calls
    
    if(write.file) {
      write.converted(hap.calls, file.name, "ab", ".txt", sep="\t")
    }
  }

  GAPIT.F = function(...) {
    # TODO http://www.zzlab.net/GAPIT/gapit_help_document.pdf
    # Required columns: 11 including rs (snp name), chrom, pos
    # Genotypes in double bit or standard iupac codes
    
    if(write.file) {
      write.converted(hap.gap, file.name, "GAPIT", ".txt", sep="\t")
    }
  }

  JOINMAP.F = function(...) {
    # TODO https://www.kyazma.nl/docs/JM4manual.pdf

    if(missing(jm.pop)) {
      stop("Joinmap population type must be specified for Joinmap export.")
    }

    check.parents(parents,hap)

    header = c("locus","segregation","phase","classification")

    # If user specifies column number instead of names
    if(is.numeric(parents)) {
      p1col = parents[1]
      p2col = parents[2]

      print(paste("Using ",colnames(hap)[parents[1]]," and ",colnames(hap)[parents[2]]," as parents.",sep=""))
    } else {
      p1col = which(grepl(parents[1],colnames(hap)))
      p2col = which(grepl(parents[2],colnames(hap)))
    }

    if(jm.pop%in% c("BC1","F2","RIx","DH","DH1","DH2","HAP","HAP1","BCpxFy","IMxFy","CM")) {
      stop("This population type is not currently implemented.")
    }

    if(jm.pop%in% c("BC1","F2","RIx","DH1","DH2","HAP1","BCpxFy","IMxFy")) {
      geno.codes = c("a","b","h","c","d","-",".","u")

    }

    if(jm.pop%in% c("DH","HAP")) {
      geno.codes = c("a","b","-",".","u")

    }

    if(any(jm.pop)=="CP") {
      pop.types = c("<abxcd>","<efxeg>","<hkxhk>","<lmxll>","<nnxnp>")
      geno.codes = c("ab","cd","ef","eg","hk","lm","ll","nn","np","--")

      # Assumes UNEAK code structure
      
    }
    
    if(write.file) {
      write.converted(hap.calls, file.name, "joinmap", ".txt", sep="\t")
    }
  }
  
  DNASP.F = function(...) {
    # TODO
    ##########################
    ## Input file for DnaSP ##
    ##########################
    hapmap3 = as.matrix(hap[as.numeric(paste(hap$present))>0.80,c(10:ncol(hap))])
    hapmap3 = hapmap3[sample(1:nrow(hapmap3), 2000),]
    hapmap3lin1 = hapmap3[,colnames(hapmap3) %in% lin1$TA]
    hapmap3lin2 = hapmap3[,colnames(hapmap3) %in% lin2$TA]
    hapmap3mc = hapmap3[,colnames(hapmap3) %in% mc$TA]
    hapmap3 = cbind(hapmap3mc, hapmap3lin1, hapmap3lin2, hapmap3)
    hapmap3 = hapmap3[,!duplicated(colnames(hapmap3))]
    hapmap3 = cbind(hapmap3, hapmap3)
    hapmap3cols_A = paste(colnames(hapmap3)[1:(ncol(hapmap3)/2)], "_A", sep="")
    hapmap3cols_B = paste(colnames(hapmap3)[1:(ncol(hapmap3)/2)], "_B", sep="")
    hapmap3cols = c(hapmap3cols_A, hapmap3cols_B)
    colnames(hapmap3) = hapmap3cols
    hapmap3[1:4,1:15]
    hapmap3[hapmap3=="H"]="N"
    hapmap3[1:4,1:15]
    hapmap3 = t(hapmap3)
    hapmap3 = cbind(rep(1:(nrow(hapmap3)/2), 2), hapmap3)
    odr = order(as.numeric(hapmap3[,1]))
    hapmap3 = hapmap3[odr,]
    hapmap3 = hapmap3[,-1]
    hapmap3 = t(hapmap3)
    hapmap3excols = cbind(rsID=paste("rs", 1:nrow(hapmap3), sep=""), position=1:nrow(hapmap3))
    hapmap3 = cbind(hapmap3excols, hapmap3)
    hapmap3[1:4,1:15]
    
    if(write.file) {
      write.converted(hapmap3, file.name, "dnasp", ".phased")
    }
  }
  
  PHYLIP.F = function(...) {
    # TODO
    phylip = as.matrix(hap[as.numeric(paste(hap$present))>0.80,10:ncol(hap)])
    phylip = t(phylip[sample(1:nrow(phylip), 2000),])
    phyliplin1 = phylip[rownames(phylip) %in% lin1$TA,]
    phyliplin2 = phylip[rownames(phylip) %in% lin2$TA,]
    phylip = rbind(phyliplin1, phyliplin2, phylip)
    phylip = phylip[!duplicated(rownames(phylip)),]
    
    phynames = paste("sequen", seq(1001,2070,by=1), sep="")
    rows = cbind(rownames(phylip), phynames)
    rownames(phylip) = phynames
    phylip[phylip=="H"]="N"
    
    if(write.file) {
      write.converted(phylip, file.name, "phylip", ".phy")
    }
  }
  
  GENO.F = function(...) {

    hap01 = hap
    hap01[,data.col:ncol(hap01)]=NA
    
    ## Define allele a and allele b
    
    if(!"alleles"%in%colnames(hap)) {
      alleles = derive.alleles(hap[,data.col:ncol(hap)])
      a = alleles$allele1
      b = alleles$allele2
    } else {
      a = substring(hap$alleles,1,1)
      b = substring(hap$alleles,3,3)
    }
    
    a[hap$alleleA<hap$alleleB] = substring(hap$alleles,3,3)[hap$alleleA<hap$alleleB]
    b[hap$alleleA<hap$alleleB] = substring(hap$alleles,1,1)[hap$alleleA<hap$alleleB]
    
    sum(a == b)
    
    # Default: turn a into -1, allele b into 1, het into 0
    hap01[hap == a] = encoding[1]
    hap01[hap == b] = encoding[3]
    hap01[hap == "H"] = encoding[2]
    
    geno = as.matrix(hap01[,data.col:ncol(hap01)])
    class(geno) = "numeric"
    geno = t(geno)
    geno
    
    if(write.file) {
      write.converted(hap.calls, file.name, "geno", ".txt", sep="\t")
    }
  }

  output = list()

  if(any(format=="MEGA")){
    output$MEGA = MEGA.F()
  }

  if(any(format=="STRUCTURE")){
    output$STRUCTURE = STRUCTURE.F()
  }

  if(any(format=="FSTRUCTURE")){
    output$FSTRUCTURE = FSTRUCTURE.F()
  }

  if(any(format=="RQTL")){
    output$RQTL = RQTL.F()
  }

  if(any(format=="AB")){
    output$AB = AB.F()
  }

  if(any(format=="GAPIT")){
    output$GAPIT = GAPIT.F()
  }

  if(any(format=="JOINMAP")){
    output$JOINMAP = JOINMAP.F()
  }
  
  if(any(format=="DNASP")){
    output$DNASP = DNASP.F()
  }
  
  if(any(format=="PHYLIP")){
    output$PHYLIP = PHYLIP.F()
  }
  
  if(any(format=="GENO")){
    output$GENO = GENO.F()
  }

  names(output) = format
  cat("Done","\n")
  invisible(output)
}

check.parents <- function(x,y) {
  if(length(x) != 2) {
    stop("Exactly two parents must be specified.")
  }

  # Check for multiple occurrences of parents, and that they both exist
  if(is.character(x)) {
    if(length(grep(paste(toupper(x),collapse="|"), toupper(colnames(y)), value=TRUE))<2) {
      stop("Unable to find both parent data columns.")
    }

    if(length(grep(paste(toupper(x),collapse="|"), toupper(colnames(y)), value=TRUE))>2) {
      stop("Found too many parent data columns. Use hap.collapse to merge.")
    }
  }
}

write.converted <- function(obj, fname, method, ext, ...) {
  write.table(obj, file=paste(fname,"_",method,".",ext,sep=""), row.names=F, quote=F, ...)
}

derive.alleles <- function(hap) {
  all.alleles = apply(hap, MARGIN = 1,function(x) names(table(unlist(lapply(x,as.character)))))
  all.alleles = lapply(alleles,function(x) x[!x%in%c("H","N")]) # TODO add option for ignorable genotypes
  
  allele1 = lapply(all.alleles, `[[`, 1)
  allele1 = unlist(allele1)
  allele2 = lapply(all.alleles, `[[`, 2)
  allele2 = unlist(allele2)
  
  alleles = data.frame(allele1,allel2)
  invisible(alleles)
}