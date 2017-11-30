# Select background variants
# Boxiang Liu (bliu2@stanford.edu)
# 2017-11-15
registerDoMC(10)

#' Load the SNPsnap dataset (https://data.broadinstitute.org/mpg/snpsnap/). If the dataset does not exist on the local computer, it will be downloaded. The file size is ~1G.  
#' 
#' @param population the population of your study
#' @return the SNPsnap dataset
#' @example
#' read_snpsnap('EUR')
read_snpsnap=function(population=c('EUR','EAS','WAFR')){
    population=match.arg(population)
    package_path=path.package('vsea')
    snpsnap_fn=sprintf('%s/inst/extdata/%s/kb1000_collection.tab.gz',package_path,population)
    
    if (!file.exists(snpsnap_fn)) {
        message('INFO - downloading SNPsnap')
        url=sprintf('https://data.broadinstitute.org/mpg/snpsnap/database/%s/kb1000/kb1000_collection.tab.gz',population)
        
        if (!dir.exists(dirname(snpsnap_fn))){
            dir.create(dirname(snpsnap_fn))
        }
        
        download.file(url,snpsnap_fn,method = "auto")
    }

    snpsnap=fread(sprintf('gunzip -c %s',snpsnap_fn),select=c(1,2,3,4,5,7,25:29))
    return(snpsnap)
}

#' Select matched background variants
#' @param x A vector of SNP IDs (format: chr:pos, e.g. 12:56115585)
#' @param snpsnap The SNPsnap dataset returned by \code{read_snpsnap}
#' @param bg_size The number of background variant per SNP ID. Default = 200
#' @return A list of \code{data.table}s. Each \code{data.table} contains the 
#' background variant for a SNP ID.
#' @details For each SNP ID, \code{select_background_variants} select 
#' \code{bg_size} number of background variants matched by minor allele 
#' frequency, distance to the nearest protein coding gene, and the number 
#' of LD proxies (r2>0.7). 
#' @example
#' x=c('12:56115585','7:18068553')
#' select_background_variants(x,snpsnap,bg_size=200)
select_background_variants=function(x,snpsnap,bg_size=200){

    background_ls=foreach(i=1:length(x))%dopar%{
        
        snpid=x[i]
        message(sprintf('INFO - %s',snpid))
        if (!snpid%in%snpsnap$snpID){
            warning('WARN - ',snpid, ' is not in the SNPsnap dataset')
            return(NA)
        }
        tmp=snpsnap[snpID==snpid,list(freq_bin,dist_nearest_gene_snpsnap_protein_coding,friends_ld07)]
        freq=tmp$freq_bin
        dist=tmp$dist_nearest_gene_snpsnap_protein_coding
        ld07=tmp$friends_ld07
        
        step=0
        tmp=data.frame()
        while (nrow(tmp)<=bg_size){
            step=step+1
            ub=1+0.1*step
            lb=1-0.1*step
            tmp=snpsnap[(freq_bin<=freq+step)&
                            (freq_bin>=freq-step)&
                            (dist_nearest_gene_snpsnap_protein_coding<=ub*dist)&
                            (dist_nearest_gene_snpsnap_protein_coding>=lb*dist)&
                            (friends_ld07<=ub*ld07)&
                            (friends_ld07>=lb*ld07),]
        }
        message(sprintf('INFO - tolerance: %s',step))
        message(sprintf('INFO - %s background variants selected',nrow(tmp)))
        
        set.seed(42)
        tmp=tmp[!snpID%in%x,]

        if (nrow(tmp)==bg_size){
            background=tmp
        } else {
            background=tryCatch(tmp[sample(1:nrow(tmp),bg_size),],error=function(e){tmp})
        }
        return(background)
    }
    names(background_ls)=x
    return(background_ls)
}

#' Merge the background_ls into an index_set
#' @param background_ls A list returned by \code{select_background_variants}
#' @param snpsnap The SNPsnap dataset returned by \code{read_snpsnap}
#' @return A \code{data.table} object
#' @details The returned \code{data.table} contains both the foreground variant 
#' and matched background variants. Each background variant is marked with the 
#' respective foreground variant in the \code{foreground_variant} column.
#' @example 
#' x=c('12:56115585','7:18068553')
#' background_ls=select_background_variants(x,snpsnap,bg_size=200)
#' index_set=create_index_set(background_ls,snpsnap)
create_index_set=function(background_ls,snpsnap){
    
    index_set=foreach(foreground_variant=names(background_ls),.combine='rbind')%dopar%{
        message('INFO - ',foreground_variant)
        if (!foreground_variant%in%snpsnap$snpID){
            message('Not in SNPsnap. Skipped.')
            return(data.table())
        }
        tmp=rbind(snpsnap[snpID==foreground_variant,],background_ls[[foreground_variant]])
        tmp$foreground_variant=foreground_variant
        return(tmp)
    }
    index_set[,chr:=str_split_fixed(snpID,':',2)[,1]]
    index_set[,pos:=as.integer(str_split_fixed(snpID,':',2)[,2])]
    return(index_set)
}