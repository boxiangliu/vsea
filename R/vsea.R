# Variant Set Enrichment Analysis
# Boxiang Liu
# 2017-11-30

#' Variant Set Enrichment Analysis
#' @param snp  A vector of SNP IDs (format: chr:pos, e.g. 12:56115585)
#' @param region A \code{data.table} of genomic regions. Must have three columns 
#' chr, start and end
#' @param plink The path to plink
#' @param bg_size The number of background variant per SNP ID. Default = 200
#' @return A list with two elements. The first element contains enrichment 
#' P-value. The second element contains enrichment odds ratio and standard error
#' @example
#' snp=c('12:56115585','7:18068553')
#' data(region)
#' res=vsea(snp,region,'EUR')
vsea=function(snp,region,population,plink='/Users/boshliu/Documents/tools/plink_mac/plink',bg_size=200){
    message('INFO - reading SNPsnap')
    snpsnap=read_snpsnap(population)
    
    message('INFO - selecting matched background variants')
    background_ls=select_background_variants(snp,snpsnap,bg_size)
    index_set=create_index_set(background_ls,snpsnap)
    
    message('INFO - adding LD proxies')
    ld_set=select_LD_variants(index_set,population,plink)
    
    message('INFO - calculating enrichment')
    overlap=overlap(ld_set,region)
    pval=calc_enrichment(overlap)
    or=calc_odds_ratio(overlap)
    return(list(pval=pval,or=or))
}