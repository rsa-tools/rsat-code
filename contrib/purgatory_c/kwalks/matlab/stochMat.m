function P=stochMat(P)

n   = size(P,1);

for i=1:n
    P(i,:) = P(i,:)/sum(P(i,:));
end    
