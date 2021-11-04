clear all
tic
vert = 150;
edg = 3000;
m = 40;%蚂蚁数量
iter_max =100;%ACA最大迭代次数
NIND = 40;%种群大小
MAXGEN = 100;%GA最高代数
AvgW=0;
TotalV = 0;
for du = 1:1
%读取图的邻接矩阵
str = strcat ('C:\Users\ASUS\Desktop\毕设\数据\mvcp3数据\newtu\n',int2str(vert),'m',int2str(edg),'gn\tu', int2str(du) , '.txt') ;
A = textread(str);
n = A(1,1);
W = A(2,:);
X = A(3:n+2,:);
edge = A(1,2);
%%定义蚁群算法参数
alpha = 3;%信息重要程度因子
beta = 2;%启发函数的重要程度因子
rho = 0.1;%信息度挥发因子
Q = 1;%常系数
Eta = 1./W;%启发函数
Tau = ones(n,n);%信息素矩阵
Table = zeros(m,n);%路径记录表
iter = 1;%迭代次数初始
Route_best = zeros(iter_max,n);%各代的最佳路径
Weight_best = zeros(iter_max,1);%各代最佳路径的长度
Weight_ave = zeros(iter_max,1);%各代路径的平均长度
%%迭代寻找最优解

while iter <= iter_max
%随机产生各个蚂蚁的起始顶点
start = zeros(m,1);

Table(:,1) = start;
%构建解空间 
vertex_index = 1:n;
%逐个蚂蚁路径选择
for i = 1:m
    %记录每只蚂蚁的运行邻接元胞矩阵
% Y = cell(m,1);
%  Y{m} = X;
 X1 = X;
temp = randperm(n);
start(i) = temp(1);
Table(i,1) = start(i);
%更新邻接矩阵
X1(start(i),:) = 0;
X1(:,start(i)) = 0;
%逐个顶点路径选择
q=2;
j=1;
    while q > 1
     j = j+1;
    q =1;
    tabu = Table(i,1:(j-1));%已访问顶点的集合
    allow_index =  ~ismember(vertex_index,tabu);
    allow = vertex_index(allow_index);%待访问顶点集合
    P = allow;
    %计算顶点间转移概率
    for k = 1:length(allow)
      P(k) = Tau(tabu(end),allow(k))^alpha+Eta(allow(k))^beta;
    end    
    P = P/sum(P);
    %轮盘赌法选择下一个顶点
    Pc = cumsum(P);
    target_index = find(Pc >= rand);
    target = allow(target_index(1));
    Table(i,j) = target;
    %更新邻接矩阵
X1(target,:) = 0;
X1(:,target) = 0;
for ii = 1:n;%求q(i)
   k=0;
    for jj = 1:n;
        k = k+X1(ii,jj);
    end
        if k <= 1
            k=0;
        end
        q = q+k;
end
    end
end
%计算各个蚂蚁的总权重
Weight = zeros(m,1);
for i = 1:m
   Route = Table(i,:);
   j=1;
      while Route(j) ~= 0
      Weight(i) = Weight(i)+W(Route(j));
      j = j+1;
       end
  
end
%计算最短权重及平均权重
if iter == 1
 [min_Weight,min_index] = min(Weight);
 Weight_best(iter) = min_Weight;
 Weight_ave(iter) = mean(Weight);
 Route_best(iter,:) = Table(min_index,:);
else
    [min_Weight,min_index] = min(Weight);
     Weight_best(iter) = min(Weight_best(iter-1),min_Weight);
     Weight_ave(iter) = mean(Weight);
    if Weight_best(iter,:) == min_Weight
        Route_best(iter,:) = Table(min_index,:);    
    else
    Route_best(iter,:) = Route_best((iter-1),:);
    end    
end

%更新信息素
Delta_Tau = zeros(n,n);
%逐个蚂蚁计算
for i = 1:m
    %逐个顶点计算
    j = 1;
   while Table(i,j+1) ~= 0
       Delta_Tau(Table(i,j),Table(i,j+1)) =  Delta_Tau(Table(i,j),Table(i,j+1))+Q/Weight(i);
       j = j+1;
   end
   j = j-1;
    Delta_Tau(Table(i,j),Table(i,1)) =  Delta_Tau(Table(i,j),Table(i,1))+Q/Weight(i);
end
Tau = (1-rho)*Tau+Delta_Tau;
%迭代次数加1，清空路径记录表
iter = iter+1;
Table = zeros(m,n);
end
%%结果展示
AA = (1:n);
Q = Route_best(iter_max,:);
for i = 1:n
AA(i) = ismember(AA(i),Q);
end

 disp(AA)
 figure(1)
 plot(1:iter_max,Weight_best,'b',1:iter_max,Weight_ave,'r')
 legend('最小权重','平均权重')
 xlabel('迭代代数')
 ylabel('权重')
title('各代最小权重与平均权重的对比')
m = A(1,2);
Pc = 0.9;%交叉概率
Pm = 0.05;%变异概率
GGAP = 0.9;%代沟
N = size(X,1);%个体染色体长度
%%定义初始化种群
Chrom = rep(AA,[NIND,1]);
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
% AvgW1 = AvgW/10
% TotalV1 = TotalV/10
% time = toc/10
% data = [vert,edg,AvgW1,TotalV1,time];
% row = size(xlsread('C:\Users\ASUS\Desktop\毕设\结果\1.xlsx'))
% str1 = strcat(int2str(row(1)+1),'A');
% xlswrite('C:\Users\ASUS\Desktop\毕设\结果\1.xlsx',data,1,str1);
