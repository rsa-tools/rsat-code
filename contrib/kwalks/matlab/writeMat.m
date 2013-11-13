function writeMat(fname,M,append)
%Usage: writeMat(fname,M,append)

nL = size(M,1);
nC = size(M,2);

if nargin < 3 
    append = 0;
end    

if append == 1
    outFile = fopen(fname,'a');
    fprintf(outFile,'#\n');
else
    outFile = fopen(fname,'w');
end    

for l=1:nL
    for c=1:nC
        fprintf(outFile,'%e',M(l,c));
        if c < nC
            fprintf(outFile,' ');
        end    
    end
    fprintf(outFile,'\n');   
end
fclose(outFile);