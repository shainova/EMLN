% NOTE!!! This file accompanies the following publication and can
% only be understood by reading the details in the manuscript and its
% SI. Please cite the original publication if using this code.
% 
% Pilosof S, Porter MA, Pascual M, Kefi S.
% The multilayer nature of ecological networks.
% Nature Ecology & Evolution (2017).


% This function analyses the modular structure of the observed temporal
% multilayer network with quantitative interlayer AND intralyer edges. The
% input is the number of realizations (runs) the function should perform
% due to the heuristic nature of the general Louvain algorithm, and the
% folder to which to write the results. It writes as output a vector (Q)
% with the maximum modularity values calculated by the function and a
% matrix (S) with module affiliations of the state nodes. The number of
% columns in S corresponds to the number of realizations and the number of
% rows is of length num_layers*num_nodes. The post processing of the output
% is done in R and expalined in the R code. The number of modules is
% calculated from S during the post processing.

function []=modularity_weighted_multilayer_obs(runs,outputfolder)

% initialize
gamma=1;
Q_runs=[];
S_runs=[];

files=dir('host_parasite_abundance_weighted_layer_*.csv'); % These are the different layers of the network, produced with the R code supplied.

% Bipartite networks have to be transformed to unipartite (square matrix).
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

N=length(A{1}); % Number of nodes (hosts+parasites)
T=length(A); % Number of layers

for r=1:runs
    r
    B=spalloc(N*T,N*T,N*N*T+2*N*T); % create an empty modularity matrix
    twomu=0;
    for s=1:T
    %%%%%% BIPARTITE:
        % In case of unipartite undirected networks the probability P of an edge exisiting between two nodes 
        % is proportional to the product of node degrees. This is
        % k_is*k_js/2m_s in eq. 1 in the SI of the paper. Note the 2*m_s because it is an
        % undirected unipartite network.
        % In the bipartite case, P is k_is*d_js/m_s. and the division is over
        % the number of edges m (see Barber 2007, eq. 15 for reasoning).
        % When networks are weighted this explanation refers to the strength
        % instead of degrees.

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
        % Here again we have to create a smymetric adjacency matrix out of the bipartite.
        onemode=zeros(p+q,p+q);
        onemode(1:p, (p+1):(p+q))=P_ij;
        onemode((p+1):(p+q), 1:p)=P_ij';
        P_ij=onemode;
        B(indx,indx)=A{s}-gamma*P_ij/m; % Note the P_ij instead of k*k' as in the unipartite version. also the m instead of 2m
    end
     
    % Because all interlayer edges have different values, OMEGA IS A MATRIX.
    % This matrix was produced using the accompanying R code.
    omega=importdata('interlayer_relative_abundance_matrix.csv');
    twomu=twomu+2*sum(sum(omega)); 
    B = B+omega;    
    
    [S,Q] = genlouvain(B,10000,0,1,1);
    Q_runs = [Q_runs Q/twomu];
    S_runs = [S_runs S];
end
dlmwrite([outputfolder,'/Q_obs.csv'],full(Q_runs)',',');
dlmwrite([outputfolder,'/S_obs.csv'],S_runs,',');

end
