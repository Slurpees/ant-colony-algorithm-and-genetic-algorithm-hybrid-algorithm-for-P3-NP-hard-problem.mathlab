clear
clc
tic
%��ȡͼ���ڽӾ���
AvgW=0;
TotalV = 0;
vert = 150;
edg = 3000;
NIND = 40;%��Ⱥ��С
MAXGEN = 200;%��ߴ���
for du = 1:1
str = strcat ('C:\Users\ASUS\Desktop\����\����\mvcp3����\newtu\n',int2str(vert),'m',int2str(edg),'gn\tu', int2str(du) , '.txt') ;
A = textread(str);
n = A(1,1);
m = A(1,2);
W = A(2,:);
X = A(3:n+2,:);
%%�����Ŵ��㷨����
Pc = 0.9;%�������
Pm = 0.05;%�������
GGAP = 0.9;%����
N = size(X,1);%����Ⱦɫ�峤��
%%�����ʼ����Ⱥ
Chrom = crtbp(NIND,N);
%%�Ż�
gen = 0;%����������
ObjV = WeightMeasure(W,Chrom);
FitnV_record = zeros(MAXGEN,1);
FitnV = Fitness(ObjV,Chrom,X);
FitnV_record(1) = mean(FitnV);
preFitnV = max(FitnV);
gen = gen+1;
while gen<MAXGEN
    %%������Ӧ��
    %fprintf('%d %1.10f\n',gen,max(FitnV))
    %line([gen-1,gen],[preFitnV,max(FitnV)]);pause(0.0001)
    %%ѡ��
    SelCh = Select(Chrom,FitnV,GGAP);
    %%�������
    SelCh = Recombin(SelCh,Pc);
    %%����
    SelCh = Mutate(SelCh,Pm);
    %%��ת����
    SelCh = Reverse(SelCh,X,W);
    %%���²����Ӵ�������Ⱥ
    ObjV1 = WeightMeasure(W,SelCh);
    FitnV1 = Fitness(ObjV1,SelCh,X);
  Chrom = Reins(Chrom,SelCh,FitnV,FitnV1);
  ObjV = WeightMeasure(W,Chrom);
  FitnV = Fitness(ObjV,Chrom,X);
  FitnV_record(gen+1,:) = mean(FitnV);
  
  %%���µ�������

  gen = gen+1;
end
%���չʾ
figure;
hold on; box on;
xlim([0,MAXGEN])
title('�Ż�����')
xlabel('����')
ylabel('����ֵ')
plot(1:MAXGEN,FitnV_record,'r')
[FitnV,index] = sort(FitnV,'descend');
AvgW = AvgW+min(ObjV);
TotalV = TotalV+sum(Chrom(index(1),:)==1);

end
AvgW1 = AvgW/10;
TotalV1 = TotalV/10;
time = toc/10;
data = [vert,edg,AvgW1,TotalV1,time];
row = size(xlsread('C:\Users\ASUS\Desktop\����\���\1.xlsx'))
str1 = strcat(int2str(row(1)+1),'A');
xlswrite('C:\Users\ASUS\Desktop\����\���\1.xlsx',data,1,str1);
