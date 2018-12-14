function writeVec(fname,N,append)
%Usage: writeN(fname,N)

n = length(N);

if nargin < 3 
    append = 0;
end    

if append == 1
    outFile = fopen(fname,'a');
    fprintf(outFile,'#\n');
else
    outFile = fopen(fname,'w');
end    

for i=1:n
    fprintf(outFile,'%e\n',N(i)); 
end
fclose(outFile);