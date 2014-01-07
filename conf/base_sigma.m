%% Conf parameters

clc;

%%%%%%%%%%%%%

% GLOBAL Conf

% always av1
% attr  _ noise _ seedType _ update _  truth _ parameter sweep _ nAgents _ forceOnV _ size 
simName = 'attrLinear_navnp_RClean_n100_fv0_s1_epsilon_v';
dumpDir = '/cluster/work/scr2/balistef/';

% we have two because we can save the new configuration in a separate
% folder analyze an old one without deleting its conf files.
CONF_SUBDIR = 'NEW2/';
SIM_SUBDIR = 'NEW2/'; 
DUMPDIR = [dumpDir SIM_SUBDIR];

%dumpDir = 'dump/';
%bsubWD = '/cluster/home/gess/balistef/matlab/hsnew/';

VIDEO = 0;
DEBUG = 0;
DUMP = 1;
DUMP_RATE = 1; % Dump every x steps
COMPUTATION = 2; % 0-local, 1-parallel, 2-LSF
SHOW_POTENTIAL = 0;
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
alphas = [0.01:0.02:0.99]; % weighting of velocity terms
Rs     = [0.01:0.01:0.1];       	   % cut-off radius
         
% ATTRACTIVE AND REPULSIVE FORCES

ks     = [1];           % Power of distance in force term

As     = [0];           % Constant in attractive force term
d0s    = [1];       	% Express the range of the interaction force (exponent divisor)

Bs     = [0];           % Constant in repulsive force term
d1s    = [1];       	% Express the range of the interaction force (exponent divisor)


% HOW EASY IS TO FIND THE TRUTH (
taus   = [1];     		% coupling coefficient (divisor)

% INDIVIDUALIZATION NOISE (position)
epsilons = 0.1;   % Std. deviation of white noise term

% NOISE on APPROACH (direction)
sigmas = [0:0.01:0.05];   % Std. deviation of white noise term

% INITIIAL VELOCITIES OF SCIENTISTS
vScalings = [0.2, 0.5, 1, 2, 10, 100]; % Scaling factor for initial (random) velocities

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

truths = [0.5; 0.5];


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

attrtype = 2;

% PLOT TYPE
plot_cross = 0;
plot_number = 1;
plot_number_color = 2;
plot_arrow = 3;

plottype = plot_cross;

% SEED TYPE
seed_fixed = 0;
seed_random = 1;
seed_seed_machinetime = 2;

seedtype = seed_random;

% NOISE TYPES
noise_on_p = 0;
noise_on_v = 1;
noise_adaptive_on_v = 2;
noise_on_v_angular = 3;
noise_navp = 4; % both 0 and 3

noisetype = 4;

% FORCES INTEGRATION on V
forces_on_v = 0;

% Seed

% This is either the fixed seed, or the seed used to init the random
% generator of random seeds before the loop of param vectorization
seed = randi(1000000);

% PARAMS 4 ANALYSIS

DUMP_ANALYSIS = 1;
DUMP_RATE_ANALYSIS = 100; % Dump every x steps

% Only clusters of size above the cutoff are included in the analysis
CLU_CUTOFF = 2;

% When computing the coverage we build a grid on top of the space of cell
% size = PRECISION
PRECISION = 100;

% Different radius to try out (around truth).
RADIUSs = [0.01, 0.05, 0.1, 0.25, 0.4];

% An agent must stay for at least STAY_FOR iterations.
STAY_FOR = 1;

% A share of at least CONSENSUS_THRESHOLD agents must stay within the
% smallest radius to be say that there is a consensus on truth.
CONSENSUS_THRESHOLD = 2/3;

% Consenus must on truth must hold for at least CONSENSUS_ON_TRUTH_FOR
% iterations to say that there is really a consensus on truth.
CONSENSUS_ON_TRUTH_FOR = 20;

% Split by Sigma

nCombinations = size(dts,2)*size(n_agents,2)*size(ideas_space_sizes,2)*...
                size(ideas_space_dims,2)*size(As,2)*size(Bs,2)*size(ks,2)*...
                size(d0s,2)*size(d1s,2)*size(alphas,2)*size(taus,2)*size(Rs,2)*...
                size(vScalings,2)*size(nClusters,2)*...
                size(clusterTightness,2)*size(truths,2)*size(attrtype,2)*...
                size(noisetype,2)*size(forces_on_v,2)*size(epsilons,2);
            
                
fprintf('%u levels of Sigma\n',  size(sigmas,2));           
fprintf('Total number of simulations = %u x %u: = %u\n', nRuns, nCombinations, nRuns*nCombinations);


% Printing parameters (has some assumptions)
fprintf('\n%s\n', simName);
fprintf('------------------------------------\n');
fprintf('%+15s = %s\n','alpha', mat2str(alphas));
fprintf('%+15s = %s\n','R', mat2str(Rs));
fprintf('%+15s = %s\n','sigma', mat2str(sigmas));
fprintf('%+15s = %s\n','epsilon', mat2str(epsilons));
fprintf('%+15s = %s\n','steps', mat2str(t_ends));
fprintf('%+15s = %s\n','nAgents', mat2str(n_agents));  
fprintf('%+15s = %s\n','IdeasSpace size', mat2str(ideas_space_sizes));
fprintf('%+15s = %s\n','tau', mat2str(taus));
fprintf('%+15s = %s\n','v_scaling', mat2str(vScalings));
fprintf('%+15s = [%2.3f:%2.3f]\n','truth', truths(1,1), truths(2,1));
fprintf('%+15s = %d\n', 'Attr. type', attrtype);
fprintf('%+15s = %d\n', 'Noise type', noisetype);
fprintf('%+15s = %d\n', 'Forces on V', forces_on_v);
fprintf('%+15s = %d\n','Seed', seed);
fprintf('------------------------------------\n');

% CREATING FOLDERS

DIR = [CONF_SUBDIR simName '/'];
if (exist(DIR, 'dir')~=0 )
    error(sprintf('Conf dir already exists: %s', DIR));
end
mkdir(DIR);

launcherMain = '../GO_FUN';
fidMain = fopen(launcherMain, 'w');

launcherCl = '../GO_ANAL_FUN';
fidCl = fopen(launcherCl, 'w');

launcherMerge = '../GO_MERGE_FUN';
fidMerge = fopen(launcherMerge, 'w');

launcherCleanup = '../GO_CLEANUP_FUN';
fidClean = fopen(launcherCleanup, 'w');

% Random seed must be initialized for each batch (level of sigma)
if (seedtype ~= seed_fixed)
    s = RandStream('mcg16807','Seed', seed);
    RandStream.setGlobalStream(s);
else
    batchSeed = seed;
end

% Save a file with all the parameters. 
% Later it will be saved one for each sigma.
paramsAllFile = sprintf('%s/params_all', DIR);
save(paramsAllFile);

old_sigmas = sigmas;
for i=1:size(sigmas,2)
    
    % Random seed must be initialized for each batch (level of sigma)
    if (seedtype ~= seed_fixed)
        batchSeed = randi(1000000);
    end 
    
    sigmas = old_sigmas(i);
    % Sigma string
    S = round(sigmas * 100);
    confFile = sprintf('%s_s%u', simName, S);
    fullName = sprintf('%s/%s', DIR, confFile);
    
    % Saving the configuration to .mat
    save(fullName);
    
    % Creating the GO_FUN file
    if (i == 1)
        cmdStr = sprintf('bsub -J hs_chain -N matlab -nodisplay -singleCompThread -r "main_fun(''conf/'',''%s'',''%s'',''%d'')"', DIR, confFile, batchSeed);
    else
        cmdStr = sprintf('bsub -J hs_chain -w "done(hs_chain)" -N matlab -nodisplay -singleCompThread -r "main_fun(''conf/'',''%s'',''%s'', ''%d'')"', DIR, confFile, batchSeed);
    end
    fprintf(fidMain, '%s\n', cmdStr);   

end

fclose(fidMain);
cmdStr = sprintf('# SCR: %s', dumpDir);
fprintf(fidCl, '%s\n', cmdStr);
cmdStr = sprintf('bsub  -R "rusage[mem=4000]" -J hs_analysis -W 24:00 -N matlab -nodisplay -singleCompThread -r "LSF_analysis(''conf/%s'')"', DIR);
fprintf(fidCl, '%s\n', cmdStr);
fclose(fidCl);

cmdStr = sprintf('bsub  -R "rusage[mem=4000]" -J hs_merge -W 24:00 -N matlab -nodisplay -singleCompThread -r "aggregate_fun(''conf/%s'')"', DIR);
fprintf(fidMerge, '%s\n', cmdStr);
fclose(fidMerge);

cmdStr = sprintf('matlab -nodisplay -singleCompThread -r "cleanup_fun(''%s%s'')"', dumpDir, DIR);
fprintf(fidClean, '%s\n', cmdStr);
fclose(fidClean);


% Saving a copy of all files in the conf dir
COPYDIR = [DIR 'launchers/']; 
mkdir(COPYDIR);
copyfile(launcherMain, COPYDIR);
copyfile(launcherCl, COPYDIR);
copyfile(launcherMerge, COPYDIR);
copyfile(launcherCleanup, COPYDIR);
