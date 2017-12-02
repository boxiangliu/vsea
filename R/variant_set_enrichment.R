# Perform VSEA
# Boxiang Liu (bliu2@stanford.edu)
# 2017-11-15

library(sinib)
#' Find overlap between variants and genomic region
#' @param ld_set The \code{data.table} returned by \code{select_LD_variants}
#' @param region A \code{data.table} of genomic regions. Must have three columns 
#' chr, start and end
#' @return A \code{data.table} that contains overlap between variants and 
#' genomic regions 
#' @example
#' set.seed(42)
#' region=make_regions(index_set,prob_fg=0.5,prob_bg=0.1)
#' overlap=overlap(ld_set,region)
overlap=function(ld_set,region){
    message('INFO - performing overlap')
    ld_set[,c('start','end'):=list(pos,pos)]
    
    stopifnot(all(c('chr','start','end')%in%names(region)))
    setkey(region,chr,start,end)
    
    # Overlap: 
    overlap=unique(foverlaps(ld_set,region[,list(chr,start,end)]))
    overlap[,c('i.start','i.end'):=NULL]
    
    overlap[,snp_overlap:=!is.na(start)]
    overlap[,loci_overlap:=any(snp_overlap),by='background_variant']
    overlap[,c('start','end'):=NULL]
    overlap=unique(overlap)
    stopifnot(nrow(overlap)==nrow(ld_set))
    stopifnot(overlap$snpID==ld_set$snpID)
    
    
    overlap=overlap[ld_proxy==FALSE,]
    overlap[,p:=mean(loci_overlap),by='foreground_variant']
    return(overlap)
}

#' Calculates enrichment p-values based on GREGOR
#' @param overlap The \code{data.table} returned by \code{overlap}
#' @return A enrichment p-value 
#' @example
#' calc_enrichment(overlap)
calc_enrichment=function(overlap){
    # Calculate enrichment p-value:
    message('INFO - calculating enrichment p-value')
    p=overlap[,list(p=unique(p)),by='foreground_variant']
    n=as.integer(rep(1,length(p$p)))
    s=sum(overlap[snpID==foreground_variant,loci_overlap])
    psinib(s,n,p$p,lower.tail=FALSE)
}

#' Calculates odds ratio between foreground variant overlap and background 
#' variant overlap (with the genomic region)
#' @param overlap The \code{data.table} returned by overlap
#' @return A \code{data.table} with mean and se of the odds ratio
#' @example
#' calc_odds_ratio(overlap)
calc_odds_ratio=function(overlap){
    message('INFO - calculating odds ratio')
    foreground_odds=overlap[background_variant==foreground_variant,sum(loci_overlap)/.N]
    background_overlap=overlap[background_variant!=foreground_variant]
    
    odds_ratio=foreach(i=1:500,.combine=c)%dopar%{
        bootstrap=background_overlap[sample(1:nrow(background_overlap),nrow(background_overlap),replace=TRUE)]
        background_odds=bootstrap[,sum(loci_overlap)/.N]
        foreground_odds/background_odds
    }
    
    mean=mean(odds_ratio)
    sd=sd(odds_ratio)
    return(data.table(mean=mean,sd=sd))
}


#' Count overlap between variant and region
#' @param overlap The object returned by overlap
#' @return The number of overlap
#' @example
#' count_overlap(overlap)
count_overlap=function(overlap){
    sum(overlap[foreground_variant==background_variant]$loci_overlap)
}