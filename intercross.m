function [a,b] = intercross(a,b)
%输入
%a,b两个待交叉的个体
%输出
%a和b为交叉后得到的两个个体
L = length(a);
r1 = randsrc(1,1,[1:L]);
r2 = randsrc(1,1,[1:L]);
    a0 = a;b0 = b;
    s = min([r1,r2]);
    e = max([r1,r2]);
    for i =s:e;
        a(i) = b0(i);
        b(i) = a0(i);
    end
    
        
