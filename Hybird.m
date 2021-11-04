clear all
tic
vert = 150;
edg = 3000;
m = 40;%��������
iter_max =100;%ACA����������
NIND = 40;%��Ⱥ��С
MAXGEN = 100;%GA��ߴ���
AvgW=0;
TotalV = 0;
for du = 1:1
%��ȡͼ���ڽӾ���
str = strcat ('C:\Users\ASUS\Desktop\����\����\mvcp3����\newtu\n',int2str(vert),'m',int2str(edg),'gn\tu', int2str(du) , '.txt') ;
A = textread(str);
n = A(1,1);
W = A(2,:);
X = A(3:n+2,:);
edge = A(1,2);
%%������Ⱥ�㷨����
alpha = 3;%��Ϣ��Ҫ�̶�����
beta = 2;%������������Ҫ�̶�����
rho = 0.1;%��Ϣ�Ȼӷ�����
Q = 1;%��ϵ��
Eta = 1./W;%��������
Tau = ones(n,n);%��Ϣ�ؾ���
Table = zeros(m,n);%·����¼��
iter = 1;%����������ʼ
Route_best = zeros(iter_max,n);%���������·��
Weight_best = zeros(iter_max,1);%�������·���ĳ���
Weight_ave = zeros(iter_max,1);%����·����ƽ������
%%����Ѱ�����Ž�

while iter <= iter_max
%��������������ϵ���ʼ����
start = zeros(m,1);

Table(:,1) = start;
%������ռ� 
vertex_index = 1:n;
%�������·��ѡ��
for i = 1:m
    %��¼ÿֻ���ϵ������ڽ�Ԫ������
% Y = cell(m,1);
%  Y{m} = X;
 X1 = X;
temp = randperm(n);
start(i) = temp(1);
Table(i,1) = start(i);
%�����ڽӾ���
X1(start(i),:) = 0;
X1(:,start(i)) = 0;
%�������·��ѡ��
q=2;
j=1;
    while q > 1
     j = j+1;
    q =1;
    tabu = Table(i,1:(j-1));%�ѷ��ʶ���ļ���
    allow_index =  ~ismember(vertex_index,tabu);
    allow = vertex_index(allow_index);%�����ʶ��㼯��
    P = allow;
    %���㶥���ת�Ƹ���
    for k = 1:length(allow)
      P(k) = Tau(tabu(end),allow(k))^alpha+Eta(allow(k))^beta;
    end    
    P = P/sum(P);
    %���̶ķ�ѡ����һ������
    Pc = cumsum(P);
    target_index = find(Pc >= rand);
    target = allow(target_index(1));
    Table(i,j) = target;
    %�����ڽӾ���
X1(target,:) = 0;
X1(:,target) = 0;
for ii = 1:n;%��q(i)
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
%����������ϵ���Ȩ��
Weight = zeros(m,1);
for i = 1:m
   Route = Table(i,:);
   j=1;
      while Route(j) ~= 0
      Weight(i) = Weight(i)+W(Route(j));
      j = j+1;
       end
  
end
%�������Ȩ�ؼ�ƽ��Ȩ��
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

%������Ϣ��
Delta_Tau = zeros(n,n);
%������ϼ���
for i = 1:m
    %����������
    j = 1;
   while Table(i,j+1) ~= 0
       Delta_Tau(Table(i,j),Table(i,j+1)) =  Delta_Tau(Table(i,j),Table(i,j+1))+Q/Weight(i);
       j = j+1;
   end
   j = j-1;
    Delta_Tau(Table(i,j),Table(i,1)) =  Delta_Tau(Table(i,j),Table(i,1))+Q/Weight(i);
end
Tau = (1-rho)*Tau+Delta_Tau;
%����������1�����·����¼��
iter = iter+1;
Table = zeros(m,n);
end
%%���չʾ
AA = (1:n);
Q = Route_best(iter_max,:);
for i = 1:n
AA(i) = ismember(AA(i),Q);
end

 disp(AA)
 figure(1)
 plot(1:iter_max,Weight_best,'b',1:iter_max,Weight_ave,'r')
 legend('��СȨ��','ƽ��Ȩ��')
 xlabel('��������')
 ylabel('Ȩ��')
title('������СȨ����ƽ��Ȩ�صĶԱ�')
m = A(1,2);
Pc = 0.9;%�������
Pm = 0.05;%�������
GGAP = 0.9;%����
N = size(X,1);%����Ⱦɫ�峤��
%%�����ʼ����Ⱥ
Chrom = rep(AA,[NIND,1]);
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
% AvgW1 = AvgW/10
% TotalV1 = TotalV/10
% time = toc/10
% data = [vert,edg,AvgW1,TotalV1,time];
% row = size(xlsread('C:\Users\ASUS\Desktop\����\���\1.xlsx'))
% str1 = strcat(int2str(row(1)+1),'A');
% xlswrite('C:\Users\ASUS\Desktop\����\���\1.xlsx',data,1,str1);
