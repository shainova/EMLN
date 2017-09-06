% NOTE!!! This file accompanies the following publication and can
% only be understood by reading the details in the manuscript and its
% SI. Please cite the original publication if using this code.
% 
% Pilosof S, Porter MA, Pascual M, Kefi S.
% The multilayer nature of ecological networks.
% Nature Ecology & Evolution (2017).


% This function is for the second null model from the publication
% (reshuffling node identities)

% The function takes as input if to permute the order of rows (=hosts) or
% columns(=parasites) by entering 1 or 2, respectively. It also receives
% the number of realizations to perform and the output folder to where
% results are written. It writes as output a vector (Q) with the maximum
% modularity values calculated by the function and a matrix (S) with module
% affiliations of the state nodes. The number of columns in S corresponds
% to the number of realizations (reshuffled networks) and the number of
% rows is of length num_layers*num_nodes. The post processing of the output
% is done in R and expalined in the R code. The number of modules is
% calculated from S during the post processing.

function []=modularity_weighted_multilayer_null2(hosts_or_parasites,runs,outputfolder)
%% initialize
gamma=1;
Q_runs=[];
S_runs=[];

files=dir('host_parasite_abundance_weighted_layer_*.csv');
omega=importdata('interlayer_relative_abundance_matrix.csv');

%% Run main loop
for r=1:runs
	r
    %% Load layers and permute
    A=cell(1,length(files));
    for i = 1:length(files)
        bip_data=importdata(files(i).name);
        [p,q]=size(bip_data);
        %% Reshuffle the node order
        % The "nodal" null model (sensu Bassett et al., 2011 PNAS)
        % reshuffles the interlayer edges between nodes and their
        % counterparts in two consecutive layers by permuting the node
        % labels separately in each layer so that node identity is not
        % preserved. In the bipartite networks, this is akin to permuting
        % the order of rows (or columns) for every given layer. You can
        % permute by rows for hosts or by columsns for parasites by
        % commenting/uncommenting the following lines:
        
        if hosts_or_parasites==1
            rowPermutations=randperm(p);
            bip_data=bip_data(rowPermutations,:);
        end
        if hosts_or_parasites==2
            colPermutations=randperm(q);
            bip_data=bip_data(:,colPermutations);
        end
        
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
    %% Create the modularity matrix
    B=spalloc(N*T,N*T,N*N*T+2*N*T);
    twomu=0;
    for s=1:T
        k=sum(A{s}); % this is the sum of degrees of the nodes in the two sets
        k=k(1:p); % This is just the first set
        d=sum(A{s}); % this is the sum of degrees of the nodes in the two sets
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
        B(indx,indx)=A{s}-gamma*P_ij/m; % Note the P_ij instead of k*k' as in the unipartite version. also the m in stead of 2m
    end
    twomu=twomu+2*sum(sum(omega)); 
    %% Run modularity
    [S,Q] = genlouvain(B,10000,0,1,1);
    Q_runs = [Q_runs Q/twomu];
    S_runs = [S_runs S];
end
%% Write results
if hosts_or_parasites==1
    dlmwrite([outputfolder,'/Q_null2_hosts.csv'],full(Q_runs)',',');
    dlmwrite([outputfolder,'/S_null2_hosts.csv'],S_runs,',');
end
if hosts_or_parasites==2
    dlmwrite([outputfolder,'/Q_null2_parasites.csv'],full(Q_runs)',',');
    dlmwrite([outputfolder,'/S_null2_parasites.csv'],S_runs,',');
end
 
end
