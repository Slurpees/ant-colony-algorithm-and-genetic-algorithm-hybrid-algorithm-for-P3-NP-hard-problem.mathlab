clear all 
clc
n = 10;
X = rand(n);
for i= 1:n
    X(i,i)=0;
end
b = reshape(X,[1,numel(X)]);
c = sort(b);
d = c(0.75*n*n);
X = X >= d
W = randi(120,1,10)

%%�����Ŵ��㷨����
NIND = 20;%��Ⱥ��С
MAXGEN = 20;%��ߴ���
Pc = 0.9;%�������
Pm = 0.05;%�������
GGAP = 0.9;%����
N = size(X,1);%����Ⱦɫ�峤��
%%�����ʼ����Ⱥ
Chrom = crtbp(NIND,N);
%%�Ż�
gen = 0;%����������
figure;
hold on; box on;
xlim([0,MAXGEN])
title('�Ż�����')
xlabel('����')
ylabel('����ֵ')
ObjV = WeightMeasure(W,Chrom)





    
