fuction ObjV = WeightMeasure1(W,Chrom)
%%����
%����ȨֵW�ͳ�ʼ��Chrom
%���
%�⼯��Ȩֵ
ObjV = zeros(NIND,1)
for j = 1:NIND
for i = 1:N
ObjV(j,1) = ObjV(j,1)+Chrom(j,i)*W[1,i];
end
end