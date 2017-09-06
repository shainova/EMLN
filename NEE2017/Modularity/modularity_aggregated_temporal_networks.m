% NOTE!!! This file accompanies the following publication and can
% only be understood by reading the details in the manuscript and its
% SI. Please cite the original publication if using this code.
% 
% Pilosof S, Porter MA, Pascual M, Kefi S.
% The multilayer nature of ecological networks.
% Nature Ecology & Evolution (2017).



% This analyses the modular structure of the observed aggregated temporal multilayer
% networks

% initialize
outputfolder='output'
runs=100;
gamma=1;

files=dir('aggregated_network_sum.csv')
A=cell(1,length(files));
X=importdata(files(1).name);
A{1}=X;
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
dlmwrite([outputfolder,'/Q_obs_aggregated_sum.csv'],full(Q_runs)',',');
dlmwrite([outputfolder,'/n_obs_aggregated_sum.csv'],full(n_runs)',',');
dlmwrite([outputfolder,'/S_obs_aggregated_sum.csv'],full(S_runs),',');


files=dir('aggregated_network_mean.csv')
A=cell(1,length(files));
X=importdata(files(1).name);
A{1}=X;
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
dlmwrite([outputfolder,'/Q_obs_aggregated_mean.csv'],full(Q_runs)',',');
dlmwrite([outputfolder,'/n_obs_aggregated_mean.csv'],full(n_runs)',',');
dlmwrite([outputfolder,'/S_obs_aggregated_mean.csv'],full(S_runs),',');
