% NOTE!!! This file accompanies the following publication and can
% only be understood by reading the details in the manuscript and its
% SI. Please cite the original publication if using this code.
% 
% Pilosof S, Porter MA, Pascual M, Kefi S.
% The multilayer nature of ecological networks.
% Nature Ecology & Evolution (2017).


% This analyses the modular structure of the observed temporal multilayer
% network when interlayer edges are zero. This procedure is different than
% the others becaues the input is single layers. So this is essentially an
% analysis of monolayer networks. This is analytically equivalent to
% analyzing a multilayer network in which all interlayer edges are zero but
% this is easier to implement because the code available for multilayer
% networks cannot produce several values of Q (one for each layer).


% initialize
outputfolder='output'
runs=100;
gamma=1;

% Choose what to analyze (uncomment/comment correctly):
files=dir('host_parasite_abundance_weighted_layer_*.csv'); % These are the different layers of the network, produced with the R code supplied.

% Load the single layers and remove hosts and parasites that do not occur
% in the layer (rows and columns that sum to zero).
A=cell(1,length(files));
for i = 1:length(files)
    X=importdata(files(i).name);
    % X(all(X==0,2),:)=[];
    % X(:,all(X==0,1))=[];
    A{i}=X;
end

T=length(A); % Number of layers

for layer=1:T
    Q_runs=[];
    n_runs=[];
    S_runs=[];
    for r=1:runs
        r
        % Modularity matrix for a bipartite monolayer network
        [B,mm]=single_layer_bipartite_B_matrix(A{layer},gamma);
        [S,Q] = genlouvain(B,10000,0,1,1);
        Q_runs = [Q_runs Q/(2*mm)];
        n_runs = [n_runs max(S)];
        S_runs = [S_runs S];
    end
    dlmwrite([outputfolder,'/Q_obs_zero_layer_',num2str(layer),'.csv'],full(Q_runs)',',');
    dlmwrite([outputfolder,'/n_obs_zero_layer_',num2str(layer),'.csv'],full(n_runs)',',');
    dlmwrite([outputfolder,'/S_obs_zero_layer_',num2str(layer),'.csv'],full(S_runs),',');
end
