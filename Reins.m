function Chrom = Reins(Chrom,SelCh,FitnV,FitnV1)
%%���²����Ӵ�������Ⱥ
%����
%Chrom ������Ⱥ
%SelCh �Ӵ���Ⱥ
%FitnV ��������Ӧ��
%���
%Chrom ��ϸ������Ӵ�������Ⱥ

NIND = size(Chrom,1);
NSel = size(SelCh,1);
Chrom = [Chrom;SelCh];
FitnV = [FitnV;FitnV1];
[FitnV,index] = sort(FitnV,'descend');
Chrom = (Chrom(index(1:NIND),:));