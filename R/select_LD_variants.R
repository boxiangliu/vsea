# Select LD variants
# Boxiang Liu (bliu2@stanford.edu)
# 2017-11-15


doMC::registerDoMC(10)

#' Select variants in LD with the index variants
#' @param index_set The output from \code{create_index_set}
#' @param population The population is one of AFR, EUR or EAS
#' @param plink The path to plink
#' @return A \code{data.table} that contains the original variants in index_set, 
#' as well as their LD proxies (r2>0.7).
#' @example
#' x=c('12:56115585','7:18068553')
#' background_ls=select_background_variants(x,snpsnap,bg_size=200)
#' index_set=create_index_set(background_ls,snpsnap)
#' select_LD_variants(index_set,'EUR','/Users/boshliu/Documents/tools/plink_mac/plink')
#' @details \code{select_LD_variants} selects LD proxies for variants in the 
#' index_set based on 1000 Genomes Phase 3. The function automatically downloads
#' 1000 Genomces Phase 3 VCF files for the specified population. The function 
#' also requires plink2 (https://www.cog-genomics.org/plink2) to perform r2 
#' calculation.

select_LD_variants=function(index_set,population,plink){
    
    tmp_dir=tempdir()
    
    package_path=path.package('vsea')
    vcf_dir=sprintf('%s/inst/extdata/vcf/',package_path)
    
    if (!dir.exists(tmp_dir)) {dir.create(tmp_dir,recursive=TRUE)}
    if (!dir.exists(vcf_dir)) {dir.create(vcf_dir,recursive=TRUE)}
    
    for (c in 1:22){
        vcf_basename=sprintf('%s.chr%s.vcf.gz',population,c)
        vcf_fn=sprintf('%s/%s',vcf_dir,vcf_basename)
        url='https://sandbox.zenodo.org/record/147424/files/'
        if (!file.exists(vcf_fn)){
            message('INFO - downloading ',vcf_fn)
            download.file(sprintf('%s/%s',url,vcf_basename),vcf_fn,method='wget')
            download.file(sprintf('%s/%s.tbi',url,vcf_basename),vcf_fn,method='wget')
        }
    }
    # Calculate LD r2:
    foreach(c=1:22)%dopar%{
        vcf_basename=sprintf('%s.chr%s.vcf.gz',population,c)
        vcf_fn=sprintf('%s/%s',vcf_dir,vcf_basename)
        
        ld_out_prefix=sprintf('%s/ld0.7.chr%s',tmp_dir,c)
        
        snp_list_fn=sprintf('%s/chr%s.snps',tmp_dir,c)
        message('INFO - writing SNP list to ',snp_list_fn)
        write.table(data.frame(unique(index_set[chr==c,snpID])),file=snp_list_fn,col.names=FALSE,row.names=FALSE,quote=FALSE)
        
        message('INFO - calculating LD for chr',c)
        command=sprintf('%s --vcf %s --keep-allele-order --r2 --ld-snp-list %s --ld-window-kb 1000 --ld-window-r2 0.7 --out %s',plink,vcf_fn,snp_list_fn,ld_out_prefix)
        print(command)
        system(command)
    }
    command=sprintf('cat %s/ld0.7.chr*.ld | grep -v "CHR_A" > %s/ld0.7.all_chr.ld',tmp_dir,tmp_dir)
    system(command)
    
    
    # Get Create LD snp set:
    message('INFO - creating LD snp set...')
    ld=fread(sprintf('%s/ld0.7.all_chr.ld',tmp_dir),col.names=c('CHR_A','BP_A','SNP_A','CHR_B','BP_B','SNP_B','R2'))
    
    ld_set=foreach(i=1:nrow(index_set),.combine='rbind')%dopar%{
        background_variant=index_set[i,snpID]
        foreground_variant=index_set[i,foreground_variant]
        
        ld_snp=data.table(snpID=ld[SNP_A==background_variant,SNP_B],r2=ld[SNP_A==background_variant,R2])
        ld_snp$background_variant=background_variant
        ld_snp$foreground_variant=foreground_variant
        ld_snp[,ld_proxy:=snpID!=background_variant]
        return(ld_snp)
    }
    
    ld_set[,chr:=paste0('chr',str_split_fixed(snpID,':',2)[,1])]
    ld_set[,pos:=as.integer(str_split_fixed(snpID,':',2)[,2])]
    
    tmp_files=list.files(tmp_dir,'^ld',full.names=TRUE)
    for (fn in tmp_files){
        file.remove(fn)
    }
    
    tmp_files=list.files(tmp_dir,'^chr',full.names=TRUE)
    for (fn in tmp_files){
        file.remove(fn)
    }
    
    return(ld_set)
}

