% NOTE!!! This file accompanies the following publication and can
% only be understood by reading the details in the manuscript and its
% SI. Please cite the original publication if using this code.
% 
% Pilosof S, Porter MA, Pascual M, Kefi S.
% The multilayer nature of ecological networks.
% Nature Ecology & Evolution (2017).


% This function is for the third null model from the publication
% (permuting the order of layers).

% The function takes as input the number of realizations to perform and the
% output folder to where results are written. It writes as output a vector
% (Q) with the maximum modularity values calculated by the function and a
% matrix (S) with module affiliations of the state nodes. The number of
% columns in S corresponds to the number of realizations (reshuffled
% networks) and the number of rows is of length num_layers*num_nodes. The
% post processing of the output is done in R and expalined in the R code.
% The number of modules is calculated from S during the post processing.

function []=modularity_weighted_multilayer_null3(runs,outputfolder)

%% initialize
gamma=1;
Q_runs=[];
S_runs=[];

files=dir('host_parasite_abundance_weighted_layer_*.csv');
omega=importdata('interlayer_relative_abundance_matrix.csv');

A=cell(1,length(files));
for i = 1:length(files)
    bip_data=importdata(files(i).name);
    % Transform the pxq matrix into (p+q)x(p+q)
    [p,q]=size(bip_data);
    onemode=zeros(p+q,p+q); 
    onemode(1:p, (p+1):(p+q))=bip_data;
    onemode((p+1):(p+q), 1:p)=bip_data';
    A{i}=onemode;
    if ~issymmetric(A{i})
        disp(['layer ',num2str(i),' is NOT symmetric'])
    end
end
N=length(A{1});
T=length(A);

%% Main loop
for r=1:runs
    r
    %% Here we permute the order of the layers
    layerPermutations=randperm(T);
    A_perm=cell(T,1);
    for x=1:T
        A_perm{x}=A{layerPermutations(x)};
    end
    
    %% Calculate modularity matrix
    B=spalloc(N*T,N*T,N*N*T+2*N*T);
    twomu=0;
    for s=1:T
        k=sum(A_perm{s}); % this is the sum of degrees of the nodes in the two sets
        k=k(1:p); % This is just the first set
        d=sum(A_perm{s}); % this is the sum of degrees of the nodes in the two sets
        d=d((p+1):(p+q)); % This is just the 2nd set
        m=sum(k); % Note the m instead of twom as in unipartite
        twomu=twomu+m; % Note that we add m and not 2m to the mu. In the unipartite version this is twomu=twomu+twom
        indx=[1:N]+(s-1)*N;
        % This calculates the matrix of probabilities accroding to eq. 15 in
        % Barber 2007
        P_ij=zeros(p,q);
        for i=1:p
            for j=1:q
                P_ij(i,j)=k(i)*d(j);
            end
        end
        % Here again we have to create a symmetric adjacency matrix out of the
        % bipartite.
        onemode=zeros(p+q,p+q);
        onemode(1:p, (p+1):(p+q))=P_ij;
        onemode((p+1):(p+q), 1:p)=P_ij';
        P_ij=onemode;
        B(indx,indx)=A_perm{s}-gamma*P_ij/m; % Note the P_ij instead of k*k' as in the unipartite version. also the m in stead of 2m
    end
 
    twomu=twomu+2*sum(sum(omega)); 
    B = B+omega;
    
    %% Run modularity
    [S,Q] = genlouvain(B,10000,0,1,1);
    Q_runs = [Q_runs Q/twomu];
    S_runs = [S_runs S];
end
%% Write results
dlmwrite([outputfolder,'/Q_null3.csv'],full(Q_runs)',',');
dlmwrite([outputfolder,'/S_null3.csv'],S_runs,',');

end
