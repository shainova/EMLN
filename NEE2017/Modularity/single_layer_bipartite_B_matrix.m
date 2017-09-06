% NOTE!!! This file accompanies the following publication and can
% only be understood by reading the details in the manuscript and its
% SI. Please cite the original publication if using this code.
% 
% Pilosof S, Porter MA, Pascual M, Kefi S.
% The multilayer nature of ecological networks.
% Nature Ecology & Evolution (2017).

function [B,mm]=single_layer_bipartite_B_matrix(A,gamma) 
    [m,n]=size(A);
    N=m+n;
    k=sum(A,2);
    d=sum(A,1);
    mm=sum(k);
    B1=A-gamma*k*d/mm;
    B=zeros(N,N);
    B(1:m,m+1:N)=B1;
    B(m+1:N,1:m)=B1';
end
