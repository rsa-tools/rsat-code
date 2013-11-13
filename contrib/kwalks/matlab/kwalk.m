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

function [N,E,DIF]=kwalk(graph,knodes,kprobas)
%Usage : [N,E,DIF]=kwalk(graph,knodes,kprobas)

k = length(knodes);

%If no initial probas, put a uniform distribution
if nargin < 3 
    kprobas = ones(k,1)/k;
end

%Check if the graph is undirected
undirected = (nargout == 3);

%Make the graph a Markov chain
n = size(graph,1);
P = stochMat(graph);
I = sparse(n,1);
I(knodes) = 1;

%Initialize the output matrices
N = sparse(n,1);
E = sparse(n,n);
if undirected == 1
    DIF = sparse(n,n);
end    

%Iterate on the nodes of interest
for i=1:k
   %Build transient and absorbing sets 
   absStates  = setdiff(knodes,knodes(i));
   tranStates = setdiff(1:n,absStates);
   nTran      = length(tranStates);
   PTran      = P(tranStates,tranStates);
   ITran      = I(tranStates);
   
   %Compute the mean passage times
   F  = inv(speye(nTran,nTran) - PTran); %Fundamental matrix     
   Z  = F'*ITran;                        %Mean Node Passage Times
   Ni = zeros(n,1);                      %Reshape in nx1   
   Ni(tranStates) = Z;          
   Ei = repmat(Ni,1,n).*P;               %Reshape in nxn
   
   %Update the output
   N = N + kprobas(i)*Ni;
   E = E + kprobas(i)*Ei; 
   
   %Update the DIF matrix for undirected graphs
   if undirected
       Ei  = kprobas(i)*Ei;
       DIF = DIF + abs(Ei-Ei');
   end    
end