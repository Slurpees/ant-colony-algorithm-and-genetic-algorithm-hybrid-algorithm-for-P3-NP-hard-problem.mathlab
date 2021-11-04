clear
clc
tic
%读取图的邻接矩阵
AvgW=0;
TotalV = 0;
vert = 150;
edg = 3000;
NIND = 40;%种群大小
MAXGEN = 200;%最高代数
for du = 1:1
str = strcat ('C:\Users\ASUS\Desktop\毕设\数据\mvcp3数据\newtu\n',int2str(vert),'m',int2str(edg),'gn\tu', int2str(du) , '.txt') ;
A = textread(str);
n = A(1,1);
m = A(1,2);
W = A(2,:);
X = A(3:n+2,:);
%%定义遗传算法参数
Pc = 0.9;%交叉概率
Pm = 0.05;%变异概率
GGAP = 0.9;%代沟
N = size(X,1);%个体染色体长度
%%定义初始化种群
Chrom = crtbp(NIND,N);
%%优化
gen = 0;%代数计数器
ObjV = WeightMeasure(W,Chrom);
FitnV_record = zeros(MAXGEN,1);
FitnV = Fitness(ObjV,Chrom,X);
FitnV_record(1) = mean(FitnV);
preFitnV = max(FitnV);
gen = gen+1;
while gen<MAXGEN
    %%计算适应度
    %fprintf('%d %1.10f\n',gen,max(FitnV))
    %line([gen-1,gen],[preFitnV,max(FitnV)]);pause(0.0001)
    %%选择
    SelCh = Select(Chrom,FitnV,GGAP);
    %%交叉操作
    SelCh = Recombin(SelCh,Pc);
    %%变异
    SelCh = Mutate(SelCh,Pm);
    %%逆转操作
    SelCh = Reverse(SelCh,X,W);
    %%重新插入子代的新种群
    ObjV1 = WeightMeasure(W,SelCh);
    FitnV1 = Fitness(ObjV1,SelCh,X);
  Chrom = Reins(Chrom,SelCh,FitnV,FitnV1);
  ObjV = WeightMeasure(W,Chrom);
  FitnV = Fitness(ObjV,Chrom,X);
  FitnV_record(gen+1,:) = mean(FitnV);
  
  %%更新迭代次数

  gen = gen+1;
end
%结果展示
figure;
hold on; box on;
xlim([0,MAXGEN])
title('优化过程')
xlabel('代数')
ylabel('最优值')
plot(1:MAXGEN,FitnV_record,'r')
[FitnV,index] = sort(FitnV,'descend');
AvgW = AvgW+min(ObjV);
TotalV = TotalV+sum(Chrom(index(1),:)==1);

end
AvgW1 = AvgW/10;
TotalV1 = TotalV/10;
time = toc/10;
data = [vert,edg,AvgW1,TotalV1,time];
row = size(xlsread('C:\Users\ASUS\Desktop\毕设\结果\1.xlsx'))
str1 = strcat(int2str(row(1)+1),'A');
xlswrite('C:\Users\ASUS\Desktop\毕设\结果\1.xlsx',data,1,str1);
