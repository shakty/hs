%% Stefano Balietti

%% Initialization of variables 
%Clear workspace
close all
clear
clc

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

% Default Values

simName = 'Sim';
dumpDir = 'dump/';
confDir = 'conf/';

compLOCAL = 0;
compPARALLEL = 1;
compLSF = 2;

%% Loading Conf
% load([confDir 'NEW/attrLinear_nv_Kseed_rndseq_tm_Rleft_n200/attrLinear_nv_Kseed_rndseq_tm_Rleft_n200_s0']);
% load([confDir 'NEW/attrLinear_nv_Kseed_rndseq_tc_Rleft/attrLinear_nv_Kseed_rndseq_tc_Rleft_s0']);
load([confDir 'NEW/attrLinear_nv_rndseed_rndseq_tm_Rleft_n100_fv0/attrLinear_nv_rndseed_rndseq_tm_Rleft_n100_fv0_s0'])

%% Modifying params locally
simName = 'RBAND';
dumpDir = '/opt/MATLAB_WORKSPACE/hs/test/'; 

VIDEO = 0;
DUMP = 1;
COMPUTATION = 0;
plottype = 0;
SHOW_POTENTIAL = 0;

% Duration
t_ends = 20;
nRuns = 1;

% Size
ideas_space_sizes = [1];
ideas_space_dims = [2];

% Scaling and nAgents
vScalings = [1];
n_agents = 100;

% Influence
alphas = 0.99;
Rs = 0.3;

% Noise
sigmas = 0.01;
epsilons = 0.1;
noisetype = 4;

% Truth
taus = 10;
truths = [0.5; 0.5];
attrtype = 2;
forces_on_v = 0;


% Initial positions


% nClusters = [0.1 0.9 0.1 0.9 ; 0.1 0.9 0.9 0.1];
% nClusters = [0.1 ; 0.1 ];

% Can contain:
%  - the coordinates of the centers of the clusters
%  - the number of clusters (1..n), centers placed randomly
%  - be equal to 0, no clusters, either init options considered
nClusters = [30]; 

% Clustered

clusterTightness = [0.05];

% Can be:
%  - equal to -1, centers are placed randomly
%  - equal to (0..n) centers are placed on a radius equal to that.
clustersInCircleOfRadius = 0.1; % [0.1:0.05:0.5]; %-1; % [-1 0.4];

% Bands
% Agents are placed randomly within a circular area (band) of area equal to
% bandArea. Inner circles will have a larger section, to make all circles
% of equal area.
% Set to -1 to avoid bands.
bandAreas = [0.4;0.5]; %computeRBands(0.1, 0.1, 0.5) % last band is slighlty smaller (0.09819; is right)

% Seed
seedtype = 0; % 0 = fixed
seed = randi(1000000); % 819325;
batchSeed = randi(1000000); % 819325;

%% Creating simName and Struct
simName = createSimName(simName,DUMP,dumpDir, 1);

simParamsStruct = struct( ...
                'dumpDir', dumpDir, ...
                'simName', simName, ...
                'VIDEO', VIDEO, ...
                'DEBUG', DEBUG, ...
                'DUMP', DUMP, ...
                'DUMP_RATE', DUMP_RATE, ...
                'nRuns', nRuns, ...
                'dts', dts, ...
                't_ends', t_ends, ...
                'n_agents', n_agents, ...
                'ideas_space_sizes', ideas_space_sizes, ...
                'ideas_space_dims', ideas_space_dims, ...
                'As', As, ...
                'Bs', Bs, ...
                'ks', ks, ...
                'd0s', d0s, ...
                'd1s', d1s, ...
                'alphas', alphas, ...
                'taus', taus, ...
                'Rs', Rs, ...
                'sigmas', sigmas, ...
                'epsilons', epsilons, ...
                'v_scalings', vScalings, ...
                'nof_clusters', nClusters, ...
                'clusterTightness', clusterTightness, ...                
                'clustersInCircleOfRadius', clustersInCircleOfRadius, ...
                'bandAreas', bandAreas, ...
                'truths', truths, ...
                'attrtype', attrtype, ...
                'noisetype', noisetype, ...
                'plottype', plottype, ...
                'SHOW_POTENTIAL', SHOW_POTENTIAL, ...
                'seedtype', seedtype, ...
                'seed', seed, ...
                'batchSeed', batchSeed, ...
                'forces_on_v', forces_on_v ...
            );


%% Store a copy of input params in DIR if DUMP is required
if (DUMP)
    struct2File( simParamsStruct, dumpDir, simName);
end

%% Start Vectorization of Parameters Sets

switch (COMPUTATION)
    case compLSF
         param_sets_LSF (simParamsStruct);
        
%    case compPARALLEL
%        param_sets_parallel (dumpDir,simName,VIDEO,DEBUG,DUMP,...
%            nRuns, dts,t_ends,n_agents,ideas_space_sizes,ideas_space_dims,...
%            As,Bs,ks,d0s,d1s,alphas,taus,Rs,sigmas,...
%            vScalings,nClusters,clusterTightness,truths);
        
    case compLOCAL
        param_sets_local (simParamsStruct);
end

fprintf('\n%s: execution completed correctly\n', simParamsStruct.simName);
% Exit Matlab when invoked from command line with -r option
%exit
