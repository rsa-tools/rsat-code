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

function res=isConnex(P)
%Usage: res=isConnex(P)

n       = size(P,1);        %Number of nodes in the graph
Np      = sparse(n,1);      %List of attained nodes
Tp      = sparse(n,1);      %Boolean array saying if a node is attained
Rr      = P;                %Boolean matrix of remaining edges    
Np(1)   = 1;                %The test start with node 1
Tp(1)   = 1;
Np_len  = 1;

%Get unused edges from nodes in Np 
[alphasF,alphasT] = find(Rr(Np(1:Np_len),:) > 0);

%Iterates until no edges can be used
while (Np_len < n) && (length(alphasF) > 0)
    %Scan the unused edges
    for i=1:length(alphasF)
        %Remove the edge
        Rr(Np(alphasF(i)),alphasT(i)) = 0;
        %If the node is not yet attained add it in Np
        if Tp(alphasT(i)) == 0
            Np_len = Np_len + 1;
            Np(Np_len) = alphasT(i);
            Tp(alphasT(i)) = 1;
        end    
    end
    %Get unused edges from nodes in Np 
    [alphasF,alphasT] = find(Rr(Np(1:Np_len),:) > 0);
end

%Check is all nodes have been reached
if Np_len == n
    res = 1;
else
    res = 0;
end    
