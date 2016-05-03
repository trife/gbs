#' Convert hap objects to other formats.
#'
#' Takes in a hap objects and converts it to formats required for other analysis programs and packages.
#'
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#'
#' @param hap The hap object to convert.
#' @param output The target format.
#' @param data.col The column corresponding to the first individual in the hap object.
#' @param write.file Should the converted object be written?
#' @param filename A name to use as the base for the output file.
#' @param parents Two columns or names of the parents being used to convert to RQTL and AB format.
#' @param encoding The numbers used to encode alleles (minor, het, major).
#' @param genotypes A vector of symbols to denote genotype calls in the hap object.
#' @param missing A vector of symbols used to denote missing data in the hap object.
#' @param pheno A data frmae of phenotypes to include in R/qtl output. See details.
#'
#' @details
#' @section Output:
#' AB, GENO, RQTL, GAPIT, STRUCTURE, FSTRUCTURE, DNASP, PHYLIP
#' 
#' \subsection{R/qtl}{
#' If the a pheno data frame is included, The first column should contain the same names that are int he hap object and each other column should be a phenotype.
#' }
#' 
#' @return
#'
#' @examples
#' data(dh)
#' dh = hap.collapse(dh,names=c("TIGER","DANBY"),match=0.9)
#' dh = filter.summary(dh,project="wheat",output="hap")$hap
#' dh.ab = hap.convert(dh,output="AB",data.col=13,parents=c("tiger","danby"))
#'
#' @export

hap.convert <- function(hap, output, data.col = 14, write.file=FALSE, filename, parents=NULL, encoding=c(-1,0,1), genotypes=c("A","C","G","T","H"), missing=c("N",NA), pheno) {
  
  # TODO add: PLINK (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed)
  # TODO add: CERVUS: three files (list possible parents, list of all progeny, genotype file)
  # TODO add: MEGA
  # TODO add: JOINMAP
  
  if(missing(output)) {
    stop("No export format specified.")
  }
  
  if(write.file==TRUE & missing(filename)) {
    stop("Must specify base filename if writing the output file.")
  }
  
  RQTL.F = function(...) {
    # Data checks
    if(length(unique(hap$chrom))<=1) {
      stop("Only one chromosome detected. Fix chrom and pos columns.")
    }
    
    if(!"chrom"%in%colnames(hap) | !"pos"%in%colnames(hap)) {
      stop("chrom or pos column missing.")
    }
    
    # Reorder based on chrom/pos
    hap = hap[with(hap, order(chrom, pos)), ]
    
    # Convert to AB
    hap.ab = hap
    
    # Remove parents
    hap.ab.nopa = hap.ab[,-which(names(hap.ab) %in% parents)]
    hap.ab.nopa = data.frame(hap.ab.nopa[,data.col:ncol(data.ab.nopa)])
    
    # Change structure
    rqtl.data = t(hap.ab.nopa)
    
    if(!missing(pheno)) {
      
      # Sort pheno data frame to match order
      pheno = pheno[match(rownames(rqtl.data),pheno[,1]),]
      
      if(!all(colnames(hap.ab.nopa) %in% pheno[,1]) || !all(rownames(rqtl.data) == pheno[,1])) {
        stop("Genotype and phenotype IDs don't match.")
      }
      
      phenotypes = names(pheno)[-1]
      numPheno = ncol(pheno) - 1
      
      rqtl.data = cbind(pheno[,-1],rqtl.data)
      
      rqtl.headers = rbind(c(rep("",numPheno),hap.ab$chrom),c(rep("",numPheno),hap.ab$pos))

      rqtl.out = rbind(rqtl.headers,rqtl.data)
      colnames(rqtl.out) = c(phenotypes,hap.ab$rs)
      
      rqtl.out = as.data.frame(rqtl.out,row.names=F)
      
      if(write.file) {
        write.converted(rqtl.out, file.name, "rqtl", ".csv", sep="\t", row.names=F)
      }
      
    } else {
      
      rqtl.headers = rbind(c("",hap.ab$chrom),c("",hap.ab$pos))
      rqtl.data = cbind(row.names(rqtl.data),rqtl.data)
      
      rqtl.out = rbind(rqtl.headers,rqtl.data)
      colnames(rqtl.out) = c("id",hap.ab$rs)
      
      rqtl.out = as.data.frame(rqtl.out,row.names=F)
      
      if(write.file) {
        write.converted(rqtl.out, file.name, "rqtl_gen", ".csv", sep="\t", row.names=F)
      }
    }
    
    rqtl.out
  }
  
  STRUCTURE.F = function(...) {
    strfile = hap
    
    # Remove duplicate tags
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
    
    # Writing out the strfile for STRUCTURE program in .txt format
    if(write.file) {
      write.converted(strfile, file.name, "structure", ".txt", sep="\t", row.names=T, col.names=F)
    }
    
    strfile
  }
  
  FSTRUCTURE.F = function(...) {
    # TODO not 100% that this is accurate

    hap.str = STRUCTURE.F()
    
    # Add 5 filler columns of empty data
    hap.fstr = cbind(rownames(hap.str),rep(NA,nrow(hap.str)),rep(NA,nrow(hap.str)),rep(NA,nrow(hap.str)),rep(NA,nrow(hap.str)),rep(NA,nrow(hap.str)),hap.str)
    
    if(write.file) {
      write.converted(hap.calls, file.name, "fstructure", ".txt", sep="\t", row.names=F)
    }
    
    hap.fstr
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
    
    hap.clean = hap[hap[,p1col]!="H" & hap[,p2col]!="H",]
    hap.clean = hap.clean[hap.clean[,p1col]!="N" & hap.clean[,p2col]!="N",]
    hap.clean = hap.clean[(hap.clean[,p1col]!=hap.clean[,p2col]),]
    
    print(paste("Removed ",nrow(hap)-nrow(hap.clean)," markers.",sep=""))
    
    ## Identify two alleles for each marker
    if(!"alleles"%in%colnames(hap.clean)) {
      alleles = derive.alleles(hap.clean[,data.col:ncol(hap.clean)])
      allele1 = alleles$allele1
      allele2 = alleles$allele2
    } else {
      allele1 = substring(hap.clean$alleles,1,1)
      allele2 = substring(hap.clean$alleles,3,3)
    }
    
    ## Identify which allele belongs to which parent and impute missing from the other parent
    hap.ab = hap.clean
    
    hap.ab[,p1col] = as.character(hap.ab[,p1col])
    hap.ab[,p2col] = as.character(hap.ab[,p2col])
    
    print("Imputing parental genotypes when only one present...")
    
    for(i in 1:nrow(hap.ab)){
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
    
    if(write.file) {
      write.converted(hap.calls, file.name, "ab", ".txt", sep="\t", row.names=F)
    }
    
    hap.calls
  }
  
  GAPIT.F = function(...) {
    gapit.out = cbind(rs = hap$rs, alleles = hap$alleles, chrom = hap$chrom, pos = hap$pos, strand=NA, assembly=NA, center=NA, alleleA=hap$alleleA, alleleB = hap$alleleB, het = hap$het, maf=hap$maf, hap[,data.col:ncol(hap)])
    
    if(write.file) {
      write.converted(hap.gap, file.name, "GAPIT", ".txt", sep="\t", row.names=F)
    }
    
    gapit.out
  }
  
  DNASP.F = function(...) {
    hapmap3 = as.matrix(hap[,data.col:ncol(hap)])
    
    # Duplicate calls and append _A and _B to column names
    hapmap3 = cbind(hapmap3, hapmap3)
    hapmap3cols_A = paste(colnames(hapmap3)[1:(ncol(hapmap3)/2)], "_A", sep="")
    hapmap3cols_B = paste(colnames(hapmap3)[1:(ncol(hapmap3)/2)], "_B", sep="")
    hapmap3cols = c(hapmap3cols_A, hapmap3cols_B)
    colnames(hapmap3) = hapmap3cols

    hapmap3[hapmap3=="H"]="N"

    hapmap3 = t(hapmap3)
    hapmap3 = cbind(rep(1:(nrow(hapmap3)/2), 2), hapmap3)
    
    odr = order(as.numeric(hapmap3[,1]))
    hapmap3 = hapmap3[odr,]
    hapmap3 = hapmap3[,-1]
    hapmap3 = t(hapmap3)
    hapmap3cols = cbind(rsID=paste("rs", hap$rs, sep=""), position=hap$pos)
    hapmap3 = cbind(hapmap3cols, hapmap3)
    
    if(write.file) {
      write.converted(hapmap3, file.name, "dnasp", ".hapmap3")
    }
    
    hapmap3
  }
  
  PHYLIP.F = function(...) {
    # TODO check for illegal characters in names: parentheses, brackets, colon, semicolon, commas
    
    phylip = as.matrix(hap[,data.col:ncol(hap)])
    phylip = t(phylip[1:nrow(phylip),])
    #phylip = phylip[!duplicated(rownames(phylip)),] TODO is this needed?
    
    phynames = 1:nrow(phylip)
    phynames = paste(phynames,rownames(phylip),sep="_")
    
    # All names must be 10 characters long, pad with spaces
    phynames = paste0(phynames,"          ")
    phynames = substr(phynames,1,10)
    
    rownames(phylip) = phynames
    
    phylip[phylip=="H"]="N"
    
    phylip = rbind(rep(NA,ncol(phylip)),phylip)
    
    phylip[1,1] = nrow(phylip) - 1
    phylip[1,2] = ncol(phylip)
    
    if(write.file) {
      write.converted(phylip, file.name, "phylip", ".phy", row.names=T, col.names=F, na="")
    }
    
    phylip
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
    
    # Default: turn minor into -1, allele major into 1, het into 0
    hap01[hap == a] = encoding[1]
    hap01[hap == b] = encoding[3]
    hap01[hap == "H"] = encoding[2]
    
    geno = as.matrix(hap01[,data.col:ncol(hap01)])
    class(geno) = "numeric"
    geno = t(geno)
    
    if(write.file) {
      write.converted(geno, file.name, "geno", ".txt", sep="\t")
    }
    
    geno
  }
  
  data.list = list()

  if(any(output=="STRUCTURE")){
    data.list$STRUCTURE = STRUCTURE.F()
  }
  
  if(any(output=="FSTRUCTURE")){
    data.list$FSTRUCTURE = FSTRUCTURE.F()
  }
  
  if(any(output=="RQTL")){
    data.list$RQTL = RQTL.F()
  }
  
  if(any(output=="AB")){
    data.list$AB = AB.F()
  }
  
  if(any(output=="GAPIT")){
    data.list$GAPIT = GAPIT.F()
  }
  
  if(any(output=="DNASP")){
    data.list$DNASP = DNASP.F()
  }
  
  if(any(output=="PHYLIP")){
    data.list$PHYLIP = PHYLIP.F()
  }
  
  if(any(output=="GENO")){
    data.list$GENO = GENO.F()
  }
  
  names(data.list) = output
  cat("Done","\n")
  invisible(data.list)
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
  write.table(obj, file=paste(fname,"_",method,ext,sep=""), quote=F, ...)
}

derive.alleles <- function(obj) {
  all.alleles = apply(obj[,data.col:ncol(obj)], MARGIN = 1,function(x) names(table(unlist(lapply(x,as.character)))))
  all.alleles = lapply(all.alleles,function(x) x[!x%in%c("H","N")]) # TODO add option for ignorable genotypes
  
  allele1 = lapply(all.alleles, `[[`, 1)
  allele1 = unlist(allele1)
  allele2 = lapply(all.alleles, `[[`, 2)
  allele2 = unlist(allele2)
  
  alleles = data.frame(allele1,allele2)
  invisible(alleles)
}