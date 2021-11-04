function Chrom = Reins(Chrom,SelCh,FitnV,FitnV1)
%%重新插入子代的新种群
%输入
%Chrom 父代种群
%SelCh 子代种群
%FitnV 父代的适应度
%输出
%Chrom 组合父代与子代的新种群

NIND = size(Chrom,1);
NSel = size(SelCh,1);
Chrom = [Chrom;SelCh];
FitnV = [FitnV;FitnV1];
[FitnV,index] = sort(FitnV,'descend');
Chrom = (Chrom(index(1:NIND),:));