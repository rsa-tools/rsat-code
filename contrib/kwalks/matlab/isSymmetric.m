%***********************************************************************%
% Author  :  Jerome Callut                                              %
% Contact :  Département d'Ingénierie Informatique                      %
%            Université catholique de Louvain                           %
%            Réaumur a.337.20                                           %
%            2, place Sainte Barbe                                      %
%            B-1348 Louvain-la-Neuve                                    %
%            tel: +32 10 47 91 16                                       %
%            e-mail : jerome.callut@info.ucl.ac.be                      %
% Modified : 03/12/2006                                                 %
%***********************************************************************%

function res=isSymmetric(A)
%Usage : res=isSymmetric(A)

n    = size(A,1);
diff = reshape(A - A',n^2,1);
idx  = find(diff ~= 0);
res  = length(idx) == 0;
