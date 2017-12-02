snpsnap=read_snpsnap(population = 'EUR')
n_peaks=2e5

region=snpsnap[sample(nrow(snpsnap),n_peaks),list(snpID,rsID)]
region[,chr:=paste0('chr',str_split_fixed(snpID,':',2)[,1])]
region[,pos:=as.integer(str_split_fixed(snpID,':',2)[,2])]
region[,c('start','end'):=list(pos,pos)]

# The probability of causal variant is P(X|Z)=a+bZ
# Z indicates whether the variant overlaps an open chromatin region.
a=0.0001
b=0.002

set.seed(42)
causal_in=rbinom(nrow(region),1,a+b)
causal_out=rbinom(nrow(snpsnap)-nrow(region),1,a)

snpsnap[,causal:=NA]
snpsnap[snpID%in%region$snpID,causal:=causal_in]
snpsnap[!snpID%in%region$snpID,causal:=causal_out]

snp=snpsnap[causal==TRUE,snpID]
res=vsea(snp,region,population = 'EUR')

