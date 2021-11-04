function SelCh = Select(Chrom,FitnV,GGAP)
NIND = size(Chrom,1);
NSel = max(floor(GGAP*NIND+ .5),2);
ChrIx = Sus(FitnV,NSel);
SelCh = Chrom(ChrIx,:);