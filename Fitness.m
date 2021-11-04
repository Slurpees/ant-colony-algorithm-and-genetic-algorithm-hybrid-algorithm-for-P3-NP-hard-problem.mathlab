function FitnV = Fitness(OobjV,Cchrom,X)
b = size(Cchrom,1);
a = size(Cchrom,2);
q = zeros(b,1);
alpha =5;
beta = 2;
%将解对应边置零
for z= 1:b%对Chrom个体关于Z循环
   X1 = X;
for m = 1:a%对个体的编码遍历
     
    if Cchrom(z,m) == 1
    X1(m,:) = 0;
    X1(:,m) = 0;
    end
end
for i = 1:a;%求q(z)
   k=0;
    for j = 1:a;
        k = k+X1(i,j);
    end
        if k <= 1
            k=0;
        end
        q(z,1) = q(z,1)+k;
end
FitnV(z,1) = 1./(OobjV(z,1)+1)+1./(1+(alpha*(q(z,1)).^beta));
end



        
        




