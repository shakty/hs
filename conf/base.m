%% Conf parameters

%%%%%%%%%%%%%

% GLOBAL Conf

%simName = 'thmiddle_av1_nv_seqrnd_attrK_tau.1';
simName = 'tec_np_seqrnd_av1_Rleft';
dumpDir = '/cluster/work/scr4/balistef/'; % dump

VIDEO = 0;
DEBUG = 0;
DUMP = 1;
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
sigmas = [0:0.1:0.5];       	% Std. deviation of white noise term

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

truths = [0.5; 0.5];

truths

% BOUNDARY CONDITIONS
bBounce = 0;
bStop   = 1;
bTorus  = 2;

boundaryCondition = bBounce;

save(simName);



nCombinations = size(dts,2)*size(n_agents,2)*size(ideas_space_sizes,2)*...
                size(ideas_space_dims,2)*size(As,2)*size(Bs,2)*size(ks,2)*...
                size(d0s,2)*size(d1s,2)*size(alphas,2)*size(taus,2)*size(Rs,2)*...
                size(sigmas,2)*size(vScalings,2)*size(nClusters,2)*...
                size(clusterTightness,2)*size(truths,2);
            
            
fprintf('Total number of simulations = %u x %u: = %u\n', nRuns, nCombinations, nRuns*nCombinations);
