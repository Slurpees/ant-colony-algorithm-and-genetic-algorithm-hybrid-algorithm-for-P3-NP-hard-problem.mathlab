function SelCh = Reverse(SelCh,X,W)
%%进化逆转函数
%输入
%SelCh 被选择的个体
%X 邻接矩阵 W 顶点权值
%输出
%SelCh 进化逆转后的个体
[row,col] = size(SelCh);
ObjV = WeightMeasure(W,SelCh);
FitnV = Fitness(ObjV,SelCh,X);
SelCh1 = SelCh;
for i = 1:row
    r1 = randsrc(1,1,[1:col]);
    r2 = randsrc(1,1,[1:col]);
    mininverse = min([r1,r2]);
    maxinverse = max([r1,r2]);
    for j = mininverse:maxinverse
    if SelCh1(i,j) == 1
        SelCh1(i,j) = 0;
    else 
        SelCh1(i,j) = 1;
    end
    end
end
    ObjV1 = WeightMeasure(W,SelCh1);
    FitnV1 = Fitness(ObjV1,SelCh1,X);
index = FitnV<FitnV1;
SelCh(index,:) = SelCh1(index,:);

    