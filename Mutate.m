function SelCh = Mutate(SelCh,Pm)
%%变异操作
%输入
%SelCh 被选择的个体
%Pm 变异概率
%输出
%SelCh 变异后的个体

[NSel,L] = size(SelCh);
for i = 1:NSel
    if Pm >rand 
        R = randperm(L);
        if SelCh(i,R(1)) == 0
            SelCh(i,R(1)) = 1;
        else SelCh(i,R(1)) = 0;
        end
    end
end