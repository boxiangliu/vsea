# Utility function
# Boxiang Liu
# 2017-11-30

#' Helper function to make a region \code{data.table}
#' @param index_set The \code{data.table} returned by \code{create_index_set}
#' @param prob The probability of sampling a region that contains a variant in
#' ld_set
make_regions=function(index_set,prob_fg,prob_bg,seed=42){
    set.seed(42)
    region=index_set[,list(snpID,chr,pos,foreground_variant)]
    if (!str_detect(region$chr[1],'chr')){
        region[,chr:=paste0('chr',chr)]
    }
    fg=region[snpID==foreground_variant]
    message(nrow(fg))
    bg=region[snpID!=foreground_variant]
    message(nrow(bg))
    fg=fg[sample(nrow(fg),round(prob_fg*nrow(fg)))]
    message(nrow(fg))
    bg=bg[sample(nrow(bg),round(prob_bg*nrow(bg)))]
    message(nrow(bg))
    region=rbind(fg,bg)
    region[,c('start','end'):=list(pos-50,pos+50)]
    return(region)
}