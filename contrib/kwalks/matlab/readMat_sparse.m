%***********************************************************************%
% Author  :  Jerome Callut                                              %
% Contact :  DŽpartement d'IngŽnierie Informatique                      %
%            UniversitŽ catholique de Louvain                           %
%            RŽaumur a.337.20                                           %
%            2, place Sainte Barbe                                      %
%            B-1348 Louvain-la-Neuve                                    %
%            tel: +32 10 47 91 16                                       %
%            e-mail : jerome.callut@info.ucl.ac.be                      %
% Modified : 03/12/2006                                                 %
%***********************************************************************%

function P=readMat_sparse(fname)
%Usage : P=readMat_sparse(fname)

inFile = fopen(fname,'r');

%Check if the file is sparse
line = fgetl(inFile);
line = textscan(line,'%n');
vec  = line{1};

%DENSE MATRIX
if length(vec) > 1
    P = load(fname,'-ascii');
%SPARSE MATRIX
else
    P = sparse(vec(1),vec(1));
    line = fgetl(inFile);
    %Parse the line number    
    while line ~= -1 
        line = textscan(line,'%s');        
        vec  = line{1};
        l    = str2num(vec{1});
        for i=2:length(vec)
            entry = textscan(vec{i},'%n:%n');
            P(l,entry{1}) = entry{2};            
        end                            
        line = fgetl(inFile);
    end 
end

fclose(inFile);