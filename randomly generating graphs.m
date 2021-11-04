clear all
n=10;
A = rand(n)
for i= 1:n
    A(i,i)=0;
end
b = reshape(A,[1,numel(A)]);
c = sort(b);
d = c(0.75*n*n);
A = A >= d;
B = randi(8,1,10)


    
