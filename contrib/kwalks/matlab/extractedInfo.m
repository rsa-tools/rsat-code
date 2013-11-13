%***********************************************************************%
% Author  :  Jerome Callut                                              %
% Contact :  Département d'Ingénierie Informatique                      %
%            Université catholique de Louvain                           %
%            Réaumur a.337.20                                           %
%            2, place Sainte Barbe                                      %
%            B-1348 Louvain-la-Neuve                                    %
%            tel: +32 10 47 91 16                                       %
%            e-mail : jerome.callut@info.ucl.ac.be                      %
% Modified : 11/12/2006                                                 %
%***********************************************************************%

function [infoU,infoC]=extractedInfo(E,K,nE)
%Usage : [infoU,infoC]=extractedInfo(E,K)

%If the subgraph is undirected only take its lower trianglar matrix
if isSymmetric(E)
    E = tril(E);
end
  
n     = size(E,1);
nK    = length(K);
tot   = sum(sum(E));
Evec  = reshape(E,n*n,1);
Epos  = find(Evec);
nPos  = length(Epos);
Evec  = full(Evec(Epos));
[V,P] = sort(Evec,'descend');

%Add the edges
sub_E    = sparse(n,n);
nodeMap  = zeros(n);
info     = zeros(nPos,7);
sub_n    = 0;
connex   = 0;
cidx     = 0;
add_K    = 0;

for i=1:nPos
    
    %Get the edge end points
    pos  = Epos(P(i)) - 1;
    from = floor(pos/n) + 1;
    to   = pos - floor(pos/n)*n + 1;
   
    %Check if the from node is already in sub_E
    if nodeMap(from) == 0
        sub_n = sub_n + 1;
        nodeMap(from) = sub_n;
        %Update the number of nodes of interrest in the subgraph
        if (add_K < nK) && (length(find(K==from)) > 0)
            add_K = add_K + 1;
        end    
    end    
    
    %Check if the from node is already in sub_E
    if nodeMap(to) == 0
        sub_n = sub_n + 1;
        nodeMap(to) = sub_n;
        %Update the number of nodes of interrest in the subgraph
        if (add_K < nK) && (length(find(K==to)) > 0)
            add_K = add_K + 1;
        end          
    end 
    
    %Add the edge in sub_E
    if from ~= to
        sub_E(nodeMap(from),nodeMap(to)) = 1;
        sub_E(nodeMap(to),nodeMap(from)) = 1;
    end
        
    %Check the graph connexity
    if connex == 0 && add_K == nK
        connex = isConnex(sub_E(1:sub_n,1:sub_n));
        if connex == 1
            cidx = i;
        end    
    end    
    
    %Add the line in the info table
    info(i,1) = i;
    info(i,2) = i/nE;
    info(i,3) = sub_n;
    info(i,4) = sub_n/n;
    info(i,5) = sum(V(1:i))/tot;
    info(i,6) = from;
    info(i,7) = to;
end

if cidx > 0
    infoU=info(1:cidx-1,:);
    infoC=info(cidx:size(info,1),:);
else
    infoU=info(:,:);
    infoC=[];
end    