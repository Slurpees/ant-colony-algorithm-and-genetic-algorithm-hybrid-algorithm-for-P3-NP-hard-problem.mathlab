%%����
%����ȨֵW�ͳ�ʼ��Chrom
%���
%�⼯��Ȩֵ
function ObjV = WeightMeasure(W,Chrom)
a = size(Chrom,1);
b = size(Chrom,2);
ObjV = zeros(a,1);
for j = 1:a
for i = 1:b
ObjV(j,1) = ObjV(j,1)+Chrom(j,i)*W(1,i);
end
end