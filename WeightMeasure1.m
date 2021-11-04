fuction ObjV = WeightMeasure1(W,Chrom)
%%输入
%顶点权值W和初始解Chrom
%输出
%解集的权值
ObjV = zeros(NIND,1)
for j = 1:NIND
for i = 1:N
ObjV(j,1) = ObjV(j,1)+Chrom(j,i)*W[1,i];
end
end