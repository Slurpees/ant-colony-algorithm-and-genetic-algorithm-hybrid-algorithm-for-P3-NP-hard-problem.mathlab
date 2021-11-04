function X = crtGraph(n)
X = rand(n);
for i= 1:n
    X(i,i)=0;
end
b = reshape(X,[1,numel(X)]);
c = sort(b);
d = c(0.75*n*n);
X = X >= d;
