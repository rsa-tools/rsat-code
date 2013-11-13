function writeMat_sparse(fname,M,append)
%Usage: writeMat_sparse(fname,M,append)  

if nargin < 3
   append = 0;
end

%SPARSE MATRIX
if issparse(M)
    if append == 1
        outFile = fopen(fname,'a');
        fprintf(outFile,'#\n');
    else
        outFile = fopen(fname,'w');
    end
    nL = size(M,1);
    fprintf(outFile,'%d\n',nL);
    for l=1:nL
        cols = find(M(l,:)~=0);
        if length(cols) > 0
            fprintf(outFile,'%d',l);
            for c=cols
                fprintf(outFile,' %d:%e',c,M(l,c)); 
            end
            fprintf(outFile,'\n');
        end                        
    end
    fclose(outFile);
%DENSE MATRIX    
else
    writeMat(fname,M,append);
end