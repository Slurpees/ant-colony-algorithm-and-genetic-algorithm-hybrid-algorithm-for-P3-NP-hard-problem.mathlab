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

%%定义遗传算法参数
NIND = 20;%种群大小
MAXGEN = 20;%最高代数
Pc = 0.9;%交叉概率
Pm = 0.05;%变异概率
GGAP = 0.9;%代沟
N = size(X,1);%个体染色体长度
%%定义初始化种群
Chrom = crtbp(NIND,N);
%%优化
gen = 0;%代数计数器
figure;
hold on; box on;
xlim([0,MAXGEN])
title('优化过程')
xlabel('代数')
ylabel('最优值')
ObjV = WeightMeasure(W,Chrom)





    
