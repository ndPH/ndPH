
% =========
% ndPH software, 
% v0.1 2017
% 
% Copyright (C) 2016, 2017 Ben Cassidy
% released under GPL license
% =========
% 
% 
% Installation:
% 1. 
% Software dependency: fda matlab package by J. Ramsay, http://www.psych.mcgill.ca/misc/fda/
% 
% 2. 
% To deal with some occasional numerical precision errors induced by extremely small and large gradients within the same curve, currently need to adjust one function in the fda package.
% 
% fdaM / smooth_monotone.m : 
% line 492 : 
% 
% if isempty(zmat)
% ->
% if false% isempty(zmat)
% 
% 3.
% Add fda and ndPH folders to Matlab path
% 
% 4. 
% Example comparing two groups of weighted networks
% 
% inputs: 
% distMat_* is a #Nodes x #Nodes x #Networks array of (ideally 
% lower-triangular) distance matrices, with elements 0 (no distance) to
% 1 (maximum distance)


rng(100) % repeatable randomness

% generate groups of random (null) networks to show that some terribly 
% ill-conditioned random partial correlation networks look similar and
% others can look different depending on data dimensionality
% 
% groups 1 & 3 can be paired (not particularly meaningful here for random)
% groups 2 & 3 have diffferent numbers of networks, must be unpaired
% groups 2 & 3 also have diffferent numbers of nodes
num_nodes = [20, 25, 20];
num_samples = [25, 50, 50];
num_networks = [10, 11, 10];

nw = cell(3,1); % cell to hold the networks
for grp = 1:3
    nw{grp} = zeros(num_nodes(grp), num_nodes(grp), num_networks(grp));
    for ii = 1:num_networks(grp)
        % static projected distance 
        % - see Cassidy 2016, "Network comparison with frequency domain
        % persistent homology"
        nw{grp}(:,:,ii) = sqrt( 1- partialcorr(randn(num_samples(grp),num_nodes(grp))).^2);
    end
end
%%
ndph_struct = cell(4,1);
%  fit ndph curves
for grp = 1:3
    ndph_struct{grp} = NDPH_fit(nw{grp}, false, 1);
end

%% plot raw Betti 0 curves
for grp = 1:3
    evalargs = ndph_struct{grp}.funcFit.evalarg_all;
    betti0_curve = num_nodes(grp) * ( 1 - eval_fd(evalargs,ndph_struct{grp}.funcFit.yhatfd_D0_all)) ;
    figure()
    plot(evalargs, betti0_curve   , 'b')
    title(['raw Betti 0 curve, group ' num2str(grp)])
end

%% plot ndph curves

for grp = 1:3
    evalargs = ndph_struct{grp}.funcFit.evalarg_all;
    ndph_curve = eval_fd(evalargs,ndph_struct{grp}.funcFit.yhatfd_D1_all_alt) ;
    figure()
    plot(evalargs, ndph_curve   , 'b')
    title(['ndph curve, group ' num2str(grp)])
end

%%
pairedORunpaired = 'unpaired';  % 'paired' or 'unpaired'
nperm = 1000;  % number of permutations
q = 0.01;  % threshold
DO_PLOT = true;
ndph_comparison_2_3 = NDPH_compare_AUtCtest( ndph_struct{2}, ndph_struct{3}, pairedORunpaired, nperm, q, DO_PLOT);
%%
pairedORunpaired = 'paired';  % 'paired' or 'unpaired'
nperm = 1000;  % number of permutations
q = 0.01;  % threshold
DO_PLOT = true;
ndph_comparison_1_3 = NDPH_compare_AUtCtest( ndph_struct{1}, ndph_struct{3}, pairedORunpaired, nperm, q, DO_PLOT);

