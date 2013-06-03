%% Stefano Balietti

%% Initialization of variables 
%Clear workspace
close all
clear
clc

%% Set Random Seed

%s = RandStream('mcg16807','Seed',0);
%RandStream.setGlobalStream(s)

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
load([confDir 'attrNormMiddle_nv_rndseq_av1_Rleft']);

%% Modifying params locally
VIDEO = 1;
DUMP = 0;
% alphas = [0.1 0.5 0.9];
COMPUTATION = 0;
attrtype = 2;
alphas = 0.1;
Rs = 0.07;
plottype = 3;
vScalings = 10;
n_agents = 10;
sigmas = 0.2;
%% Creating simName and Struct
simName = createSimName(simName,DUMP,dumpDir);

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
                'v_scalings', vScalings, ...
                'nof_clusters', nClusters, ...
                'clusterTightness', clusterTightness, ...
                'truths', truths, ...
                'attrtype', attrtype, ...
                'noisetype', noisetype, ...
                'plottype', plottype, ...
                'seedtype', seed_fixed ...
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
