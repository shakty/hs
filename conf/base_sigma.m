%% Conf parameters

clc;

%%%%%%%%%%%%%

% GLOBAL Conf

% always av1
% attr  _ noise _ seedType _ update _  truth _ parameter sweep
simName = 'attrFunnel_nv_Kseed_rndseq_tc_Rleft';
dumpDir = '/cluster/work/scr5/balistef/'; 
%dumpDir = 'dump/';
bsubWD = '/cluster/home/gess/balistef/matlab/hsnew/';

VIDEO = 0;
DEBUG = 0;
DUMP = 1;
DUMP_RATE = 1; % Dump every x steps
COMPUTATION = 2; % 0-local, 1-parallel, 2-LSF

%%%%%%%%%%%%%

% MODEL Conf

nRuns = 1;             % Number of simulation runs with same param set

dts = [0.01];           % time_step
t_ends = [20];          % running time

n_agents = [100];       % number of agents

ideas_space_sizes = [1];% size of ideas space
ideas_space_dims = [2]; % dimension of ideas space

% If A = B repulsion and attraction nullify
% If A > B and there is some initial speed -> speedy bounce is possible at boundaries

% alpha is between 0 and 1. 
%
% 0 -> velocity is just from the average within the radius;
% 1-> just your own

% ks the bigger the less groups

% VELOCITY 
alphas = [0:0.01:1];       	% weighting of velocity terms
Rs     = [0:0.01:0.3];       	% cut-off radius

% ATTRACTIVE AND REPULSIVE FORCES

ks     = [1];           % Power of distance in force term

As     = [0];           % Constant in attractive force term
d0s    = [1];       	% Express the range of the interaction force (exponent divisor)

Bs     = [0];           % Constant in repulsive force term
d1s    = [1];       	% Express the range of the interaction force (exponent divisor)


% HOW EASY IS TO FIND THE TRUTH (
taus   = [0.1];     		% coupling coefficient (divisor)

% WHITE NOISE
sigmas = [0:0.1:0.5];   % Std. deviation of white noise term

% INITIIAL VELOCITIES OF SCIENTISTS
vScalings = [1];     	% Scaling factor for initial (random) velocities

% INITIAL POSITIONS OF SCIENTISTS
nClusters = [0];    	% number of clusters of the initial positions
clusterTightness = [0.25]; % Tightness of clusters

% TRUTH POSITION
% Generate Truth Vector for 2D Truth

%hGrid = [0:0.25:1];
hGrid = 0.25:0.25:0.75;
nPointsGrid = length(hGrid);
vGrid = zeros(nPointsGrid,1)';
for i=2:numel(hGrid)
    vGrid = [ vGrid repmat(hGrid(i),1,nPointsGrid)];
end
truths = [repmat(hGrid,1,nPointsGrid); vGrid];

truths = [0.1; 0.1];


% BOUNDARY CONDITIONS (not used yet)
bBounce = 0;
bStop   = 1;
bTorus  = 2;
boundaryCondition = bBounce;


% TRUTH ATTRACTION FORCE TYPE
attr_zero = 0;
attr_const = 1;
attr_linear = 2;
attr_expo = 3;
attr_millean_arena = 4;
attr_hard_to_find = 5;
attr_wide_funnel = 6;
attr_gentle_landing = 7;

attrtype = 6;

% PLOT TYPE
plot_cross = 0;
plot_number = 1;
plot_number_color = 2;
plot_arrow = 3;

plottype = plot_cross;

% SEED TYPE
seed_fixed = 0;
seed_random = 1;

seedtype = 0;

% NOISE TYPES
noise_on_p = 0;
noise_on_v = 1;
noise_adaptive_on_v = 2;

noisetype = 1;


% Split by Sigma

nCombinations = size(dts,2)*size(n_agents,2)*size(ideas_space_sizes,2)*...
                size(ideas_space_dims,2)*size(As,2)*size(Bs,2)*size(ks,2)*...
                size(d0s,2)*size(d1s,2)*size(alphas,2)*size(taus,2)*size(Rs,2)*...
                size(vScalings,2)*size(nClusters,2)*...
                size(clusterTightness,2)*size(truths,2)*size(attrtype,2)*...
                size(noisetype,2);
            
                
fprintf('%u levels of Sigma\n',  size(sigmas,2));           
fprintf('Total number of simulations = %u x %u: = %u\n', nRuns, nCombinations, nRuns*nCombinations);


% CREATING FOLDERS
CONF_SUBDIR = 'NEW/';
DIR = [CONF_SUBDIR simName '/'];
if (exist(DIR, 'dir')~=0 )
    error('Dir already exists');
end
mkdir(DIR);

launcherMain = '../GO_FUN';
fidMain = fopen(launcherMain, 'w');

launcherCl = '../GOCL_FUN';
fidCl = fopen(launcherCl, 'w');


file_merge = '../bash_merge_csv';
fidFileMerge = fopen(file_merge, 'w');
cmdStr = sprintf('#!/bin/sh');
fprintf(fidFileMerge, '%s\n', cmdStr);
cmdStr = sprintf('OUTFILE_PARAMS="%s%s%s"', dumpDir, DIR, 'params_all.csv');
fprintf(fidFileMerge, '%s\n', cmdStr);
cmdStr = sprintf('OUTFILE_MACRO="%s%s%s"', dumpDir, DIR, 'clusters_macro_all.csv');
fprintf(fidFileMerge, '%s\n', cmdStr);
cmdStr = sprintf('OUTFILE_MICRO="%s%s%s"', dumpDir, DIR, 'clusters_micro_all.csv');
fprintf(fidFileMerge, '%s\n', cmdStr);
cmdStr = sprintf('OUTFILE_MACRO_AVG_SPLIT="%s%s%s"', dumpDir, DIR, 'clusters_macro_avg_split.csv');
fprintf(fidFileMerge, '%s\n', cmdStr);


old_sigmas = sigmas;
for i=1:size(sigmas,2)

    sigmas = old_sigmas(i);
    % Sigma string
    S = sigmas*10;
    confFile = sprintf('%s_s%u', simName, S);
    fullName = sprintf('%s/%s', DIR, confFile);
    
    % Saving the configuration to .mat
    save(fullName);
    
    % Creating the GO_FUN file
    if (i == 1)
        cmdStr = sprintf('bsub -J hs_chain -W 36:00 -N matlab -nodisplay -singleCompThread -r "main_fun(''conf/'',''%s'',''%s'')"', DIR, confFile);
    else
        cmdStr = sprintf('bsub -J hs_chain -w "done(hs_chain)" -W 36:00 -N matlab -nodisplay -singleCompThread -r "main_fun(''conf/'',''%s'',''%s'')"', DIR, confFile);
    end
    fprintf(fidMain, '%s\n', cmdStr);   
    
    % Creating the GOCL_FUN
    cmdStr = sprintf('bsub -J hs_cl_%u -W 36:00 -N matlab -R "rusage[mem=20000]" -nodisplay -singleCompThread -r "temporalysis_fun(''%s'',''%s'',''%s'')"', S, dumpDir, DIR, confFile);
    fprintf(fidCl, '%s\n', cmdStr);
    
    % Creating bash_merge_csv
    if (i == 1)
        cmdStr = sprintf('cat %s%s%s/%s > $OUTFILE_PARAMS', dumpDir, DIR, confFile, 'params.csv');       
        cmdStr1 = sprintf('cat %s%s%s/%s > $OUTFILE_MACRO', dumpDir, DIR, confFile, 'clusters_macro.csv');
        cmdStr2 = sprintf('cat %s%s%s/%s > $OUTFILE_MICRO', dumpDir, DIR, confFile, 'clusters_micro.csv');
        cmdStr3 = sprintf('cat %s%s%s/%s > $OUTFILE_MACRO_AVG_SPLIT', dumpDir, DIR, confFile, 'clusters_macro_avg.csv');
    else
        cmdStr = sprintf('sed -e ''1d'' %s%s%s/%s >> $OUTFILE_PARAMS', dumpDir, DIR, confFile, 'params.csv');
        cmdStr1 = sprintf('sed -e ''1d'' %s%s%s/%s >> $OUTFILE_MACRO', dumpDir, DIR, confFile, 'clusters_macro.csv');
        cmdStr2 = sprintf('sed -e ''1d'' %s%s%s/%s >> $OUTFILE_MICRO', dumpDir, DIR, confFile, 'clusters_micro.csv');
        cmdStr3 = sprintf('sed -e ''1d'' %s%s%s/%s >> $OUTFILE_MACRO_AVG_SPLIT', dumpDir, DIR, confFile, 'clusters_macro_avg.csv');
    end
    fprintf(fidFileMerge, '%s\n', cmdStr);
    fprintf(fidFileMerge, '%s\n', cmdStr1);
    fprintf(fidFileMerge, '%s\n', cmdStr2);
    fprintf(fidFileMerge, '%s\n', cmdStr3);
end
fclose(fidMain);
fclose(fidCl);
fclose(fidFileMerge);

% Creating the GOMERGECSV_FUN
launcherMerge = '../GOMERGECSV_FUN';
fidMerge = fopen(launcherMerge, 'w');
cmdStr = sprintf('bsub -J hs_csv_merge -W 1:00 -N %sbash_merge_csv', bsubWD);
fprintf(fidMerge, '%s\n', cmdStr);
fclose(fidMerge);


% Creating the GOAGGR_FUN
launcherAggr = '../GOAGGR_FUN';
fidAggr = fopen(launcherAggr, 'w');
cmdStr = sprintf('bsub -J hs_cl_aggr -W 1:00 -N matlab -nodisplay -singleCompThread -r "aggregate_fun(''%s'',''%s'')"', dumpDir, DIR );
fprintf(fidAggr, '%s\n', cmdStr);
fclose(fidAggr);

% Saving a copy of all files in the conf dir
COPYDIR = [DIR 'launchers/']; 
mkdir(COPYDIR);
copyfile(launcherMerge, COPYDIR);
copyfile(launcherAggr, COPYDIR);
copyfile(file_merge, COPYDIR);
copyfile(launcherMain, COPYDIR);
copyfile(launcherCl, COPYDIR);
