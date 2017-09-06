% NOTE!!! This file accompanies the following publication and can
% only be understood by reading the details in the manuscript and its
% SI. Please cite the original publication if using this code.
% 
% Pilosof S, Porter MA, Pascual M, Kefi S.
% The multilayer nature of ecological networks.
% Nature Ecology & Evolution (2017).


%Create modular matrices. Block1 is the plant-herbivore network and Block2
%is the plant-parasite network.
blockMatrix1=zeros(9,9);
blockMatrix1(1:3,1:3)=1;
blockMatrix1(4:6,4:6)=1;
blockMatrix1(7:9,7:9)=1;

blockMatrix2=zeros(9,9);
blockMatrix2(1:2,1:3)=1;
blockMatrix2(3:3,4:6)=1;
blockMatrix2(4:4,7:9)=1;
blockMatrix2(5:5,4:6)=1;
blockMatrix2(6:6,1:3)=1;
blockMatrix2(7:7,7:9)=1;
blockMatrix2(8:8,4:6)=1;
blockMatrix2(9:9,7:9)=1;

% Plot the networks. See the paper by Flores et al. 2015 on how to install
% BiMat which allows plotting the bipartite networks.
bp = Bipartite(blockMatrix1)
bp.plotter.use_type_interaction = true;
bp.plotter.color_interactions(1,:) = [0 0 0]; %Red color for clear lysis
bp.plotter.color_interactions(2,:) = [0 0 0]; %Blue color for turbid spots
bp.plotter.back_color = 'white';
bp1 = Bipartite(blockMatrix2)
bp1.plotter.use_type_interaction = true;
bp1.plotter.color_interactions(1,:) = [0 0 0]; %Red color for clear lysis
bp1.plotter.color_interactions(2,:) = [0 0 0]; %Blue color for turbid spots
bp1.plotter.back_color = 'white';

figure(1);
set(gcf,'Position',[0+100 72 1754 922]);
subplot(2,2,2);
bp1.plotter.PlotGraph();
subplot(2,2,4);
bp1.plotter.PlotMatrix();
subplot(2,2,1);
bp.plotter.PlotGraph();
set(gca,'xdir','reverse');
subplot(2,2,3);
bp.plotter.PlotMatrix();
set(gca,'xdir','reverse');


% Calculate modularity in each layer separately.
[B,mm]=single_layer_bipartite_B_matrix(blockMatrix1,1);
[S,Q]=genlouvain(B,10000,0,1,1);
Q/(2*mm)
S
[B,mm]=single_layer_bipartite_B_matrix(blockMatrix2,1);
[S,Q]=genlouvain(B,10000,0,1,1);
Q/(2*mm)
S

%Multilayer modularity
runs=1000;
Omegas=[0 0.5 1000]; % These values correspond to the figure in the paper
S_sweep=zeros(36, runs, length(Omegas)); % 36 is the total number of STATE nodes
Q_sweep=zeros(runs, length(Omegas));
n_sweep=zeros(runs, length(Omegas));
moduleOverlap=zeros(9, runs, length(Omegas)); % 9 is the number of plants

for w=1:length(Omegas)
    omega=Omegas(w);
    [B_multilayer,twomu] = bipartite_modularity_diag_coupling(blockMatrix1,blockMatrix2,1,omega,0);
    for i=1:runs
        [Sw,Qw] = genlouvain(B_multilayer,10000,0,1,1);
        S_sweep(:,i,w)=Sw;
        Q_sweep(i,w)=Qw/twomu;
        n_sweep(i,w)=max(Sw);
        plants1=Sw(1:9);
        plants2=Sw(19:27);
        moduleOverlap(:,i,w)=plants1-plants2;
    end
end

% calculate the number of plants for which the state nodes were assigned to
% different modules
for w=1:length(Omegas)
    mean(sum(moduleOverlap(:,:,w)==0))
end

% Find the run with maximum modularity
maxQ=find(Q_sweep(:,3)==max(Q_sweep(:,3)))
maxQ=maxQ(1)
S_sweep(:,maxQ,3)


