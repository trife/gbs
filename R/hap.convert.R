#' Convert hap objects to other formats
#' 
#' Takes in hap files and converts them to the formats needed for other analysis programs
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Traci Kantarski, \email{tkantarski@@ksu.edu}
#' 
#' @param hap.obj the hap object to convert
#' @param format the format you wish to convert the hap object to
#' @param parents columns for the two parents being used to convert to RQTL and AB format
#' 
#' @keywords 
#' 
#' @examples
#' 
#' @export

hap.convert <- function(hap,format=c("MEGA","STRUCTURE","FSTRUCTURE","RQTL","AB","GAPIT","JOINMAP"),parents=NULL,jm.pop = c("BC1","F2","RIx","DH","DH1","DH2","HAP","HAP1","CP","BCpxFy","IMxFy")) {
  
  # TODO add: PLINK (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed), CERVUS, DNASP
  
  if(missing(format)) {
    stop("No export format specified.")
  }

  MEGA.F = function(...) {
    # TODO http://www.megasoftware.net/mega4/mega4.pdf
  }
  
  RQTL.F = function(...) {
    # TODO http://www.inside-r.org/packages/cran/qtl/docs/read.cross
    
    # line id in first column, marker names across header, chromosome on second line, A B H genotypes
  }
  
  STRUCTURE.F = function(...) {
    # TODO http://pritchardlab.stanford.edu/structure_software/release_versions/v2.3.4/structure_doc.pdf
  }
  
  FSTRUCTURE.F = function(...) {
    # TODO https://github.com/rajanil/fastStructure
    
    # rows in the data file correspond to samples, with two rows per sample and columns correspond to SNPs
    # The first 6 columns of the file will be ignored; typically include IDs, metadata, etc. 
    # only handles bi-allelic loci. The two alleles at each locus can be encoded as desired
    # missing data should be encoded as -9.
  }
  
  AB.F = function(...) {
    check.parents(parents,hap)

    # If user specifies column number instead of names
    if(is.numeric(parents)) {
      p1col = parents[1]
      p2col = parents[2]
      
      print(paste("Using ",colnames(hap)[parents[1]]," and ",colnames(hap)[parents[2]]," as parents.",sep=""))
    } else {
      p1col = which(grepl(parents[1],colnames(hap)))
      p2col = which(grepl(parents[2],colnames(hap)))
    }
    
    ## Remove markers that are hets, identical, or missing in both parents
    print("Removing markers that are het, identical, or missing in both parents...")

    hap.clean = hap[hap[,p1col]!="H" & hap[,p2col]!="H",]
    hap.clean = hap.clean[hap.clean[,p1col]!="N" & hap.clean[,p2col]!="N",]
    hap.clean = hap.clean[(hap.clean[,p1col]!=hap.clean[,p2col]),]

    print(paste("Removed ",nrow(hap)-nrow(hap.clean)," markers.",sep=""))
    
    ## Identify two alleles for each marker
    alleles = apply(hap.clean,MARGIN = 1,function(x) names(table(unlist(lapply(x,as.character)))))
    alleles = lapply(alleles,function(x) x[!x%in%c("H","N")]) # TODO add option for ignorable genotypes
    
    ## Identify which allele belongs to which parent and impute missing from the other parent
    # TODO add an if here to skip this if there's no missing data
    p1 = hap.clean[,p1col]
    p2 = hap.clean[,p2col]
    
    allele1 = lapply(alleles, `[[`, 1)
    allele1 = unlist(allele1)
    allele2 = lapply(alleles, `[[`, 2)
    allele2 = unlist(allele2)
    
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
  }
  
  GAPIT.F = function(...) {
    # TODO http://www.zzlab.net/GAPIT/gapit_help_document.pdf
    # Required columns: 11 including rs (snp name), chrom, pos
    # Genotypes in double bit or standard iupac codes
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

      # Custom Genotype Calling & Data Filtering
      
      call.hets <- function(x)  {
        Counts <- as.numeric(strsplit(as.character(x),"\\|")[[1]])
        Total <- sum(Counts)
        Minor <- min(Counts)
        Major <- which(Counts==max(Counts))[1]
        R1 <- Minor/Total
        {    
          
          # If there are no reads, call it NA
          
          if (Total == 0)	{
            return(NA)
          }
          
          # This is where genotype calls are customized.
          # If one allele has no counts, the second allele must have 5 counts or more to be called homozygous for that allele.
          # If there is less than 5 counts, NA is returned 
          
          else if (Minor == 0 && Total <= 4) {
            return(NA)
          }
          
          else if (Minor == 0 && Total > 4)	{
            return(c(0,2)[Major])
          }
          
          # If the minor allele frequency for a given genotype is less than 10%,  minor allele is considered a likely error
          # and the genotype is considered homozygous for the major allele
          
          else if (Minor > 0 &&  R1 < 0.10)	{
            return(c(0,2)[Major])
          }
          
          else if (Minor > 0 &&  R1 >= 0.10)	{
            return(1)
          }
          else return(1)
        }
      }
      
      count.hets <- function(Matrix)	{
        Mat <- as.matrix(Matrix)
        Vec <- as.vector(Mat)
        
        out.Vec <- sapply(Vec, call.hets)
        
        out.Mat <- matrix(out.Vec, nrow=nrow(Mat))
        return(out.Mat)
      }
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
  
  names(output) = format
  cat("DONE","\n")
  invisible(output)
}

check.parents <- function(x,y) {
  if(length(x) != 2) {
    stop("Exactly two parents must be specified.")
  }
  
  # Check for multiple occurrences of parents, and that they both exist
  if(is.character(x)) {
    if(length(grep(paste(x,collapse="|"), colnames(y), value=TRUE))<2) {
      stop("Unable to find both parent data columns.")
    }
    
    if(length(grep(paste(x,collapse="|"), colnames(y), value=TRUE))>2) {
      stop("Found too many parent data columns. Use hap.collapse to merge.")
    }
  }
}