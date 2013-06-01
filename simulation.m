function conv = simulation (params)

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

%params

folderName = params.folderName;
simCount = params.simCount;
VIDEO = params.VIDEO;
DEBUG = params.DEBUG;
DUMP = params.DUMP;
DUMP_RATE = params.DUMP_RATE;
nRun = params.nRun;
dt = params.dt;
t_end = params.t_end;
n_agents = params.n_agents;
ideas_space_size = params.ideas_space_size;
ideas_space_dim = params.ideas_space_dim;
A = params.A;
B = params.B;
k = params.k;
d0 = params.d0;
d1 = params.d1;
alpha = params.alpha;
tau = params.tau;
R = params.R;
sigma = params.sigma;
v_scaling = params.v_scaling;
nof_cluster = params.nof_cluster;
clusterTightness = params.clusterTightness;
truth = params.truth;

display('Starting...');

% SEED TYPE
seed_fixed = 0;
seed_random = 1;

if (params.seedtype == seed_fixed)
    % Set random seed
    s = RandStream('mcg16807','Seed', 0);
    RandStream.setGlobalStream(s);
end



%% Initializations and Definitions

% Row-vector holding agents' positions.
% agents = ideas_space_size.*rand(ideas_space_dim,n_agents);
agents = initial_pos_clustered(nof_cluster,clusterTightness,n_agents,...
                                ideas_space_size,ideas_space_dim);

colors = {'magenta','yellow','black', 'cyan', 'red', 'green', 'blue'};
  
% Row-vector holding agents' velocities
v = v_scaling*(rand(ideas_space_dim,n_agents)-0.5);

rownorm = @(X,P) sum(abs(X).^P,2).^(1/P);
colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);

% PLOT TYPE
plot_cross = 0;
plot_number = 1;
plot_number_color = 1;

% Attraction to truth force type
attr_zero = 0;
attr_const = 1;
attr_linear = 2;
attr_normal_middle = 3;
attr_normal_closer_t = 4;
attr_lognormal = 5;

switch (params.attrtype)
    case attr_zero
    % NO TRUTH
    ths = @(x) (zeros(2, n_agents));
    
    case attr_const
    % TRUTH Constant
    ths = @(x) (repmat(truth,1,n_agents)-x).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-x,2),2,1));

    case attr_linear
    % TRUTH Linear
    ths = @(x) (repmat(truth,1,n_agents)-x)./tau;
    
    case attr_normal_middle
    % TRUTH Accelerating in the middle (NORMAL)
    ths = @(x) (repmat(truth,1,n_agents)-x).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-x,2),2,1)).*normpdf(abs(repmat(truth,1,n_agents)-x),norm(truth)./2,0.1);

    case attr_normal_closer_t
    % TRUTH Accelerating in the middle (NORMAL) mean closer to TRUTH
    ths = @(x) (repmat(truth,1,n_agents)-x).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-x,2),2,1)).*normpdf(abs(repmat(truth,1,n_agents)-x),norm(truth)./3,0.1);

    case attr_lognormal
    % TRUTH Accelerating closer to Truth (LOG-NORMAL)
    ths = @(x) (repmat(truth,1,n_agents)-x).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-x,2),2,1)).*lognpdf(abs(repmat(truth,1,n_agents)-x),log(norm(truth)),1);
end

% NOISE TYPES
noise_on_p = 0;
noise_on_v = 1;
noise_adaptive_on_v =2;


% d: Eucledian distance points i,j
d = @(agents,i,j) norm(agents(:,i)-agents(:,j));  

% df: Attraction force between agents i and j
df = @(agents,i,j) ( A*d(agents,i,j)^k * exp( -d(agents,i,j) / d0) - B*d(agents,i,j)^k * exp( -d(agents,i,j) / d1)) .* ( (agents(:,j) - agents(:,i)) / d(agents,j,i) );                 

%% Initial Plot
if (VIDEO)
    fig1 = figure(1);

    % Initial Positioning
    hold on
    plot(agents(1,:),agents(2,:),'rx');
    plot(truth(1),truth(2),'go');
    hold off
end


if (DUMP)
    agentspos = agents;
    agentsv = v;
end

%% Time iteration

% Start Timer
tStart = tic; 

 
for t=0:dt:t_end                      
    
    %% Compute Total Attraction Force
    
    % agents are columns (i), dims stored in rows, sum every entry over j
    force = zeros(ideas_space_dim,n_agents);                     
    if (A ~= 0 && A ~= B) % either A or B is ~- 0 and they are not both equal
    for i=1:n_agents
        % Generates the ids of all the other agents
        other_agents_ids = mod([i:n_agents+i-2],n_agents)+1;
        for ii=1:size(other_agents_ids)
            j = other_agents_ids(ii);
            force_current = df(agents,i,j);
            force(:,i)= force(:,i) + force_current;
        end
    end
    end
     
    
    %% Compute Velocity influenced by agents in cutoff
    for i=1:n_agents
        % Important: the update is sequential and random in the order
        other_agents_ids = mod([i:n_agents+i-2],n_agents)+1;
        other_agents_ids = other_agents_ids(randperm(length(other_agents_ids)));
        agents_in_cutoff=zeros(1,n_agents);
        for ii=1:size(other_agents_ids,2)
            j = other_agents_ids(ii);
            if (d(agents,i,j) <= R)
                agents_in_cutoff(j)=1;
                %v(:,i)=v(:,i)+(1-alpha).*av_velocity(i,R); 
                %agents_in_cutoff=find(d(i,j)<R);
            end
        end

        % Important: only if some agents are in the cutoff radius
        % the velocity is split between alpha and 1-alpha.
        % If not the agent continuous his search undisturbed
        if (nnz(agents_in_cutoff) ~=0) 
            average_velocity=(v*(agents_in_cutoff)') ./ nnz(agents_in_cutoff);  
            v(:,i) = alpha * v(:,i) + (1-alpha) * average_velocity;   
        end
       
        
        if (params.noisetype == noise_on_v)
            % Important: noise on velocity
            v(:,i) = v(:,i) + normrnd(0,sigma,ideas_space_dim,1);
            
        else
            if (params.noisetype == noise_adaptive_on_v)
            % Adaptive noise: the faster you go the higher the noise
            v(:,i) = normrnd(v(:,i),sigma,ideas_space_dim,1);
            end
        end
            
        
    end
    
    %% Boundary Conditions: reflecting from the walls
    for n=1:n_agents
    	for m=1:ideas_space_dim 
            %if (agents(m,n)<=0 || agents(m,n)>=ideas_space_size)
            %   v(m,n)=-v(m,n);
            %end
            if (agents(m,n) <= 0)
                %agents(m,n)=0;     %ugly but effective. Is this okay?
                
                % Important Difference!!! whether Bouncing or Zero v
                %v(m,n)=0;
                v(m,n)=abs(v(m,n)); %avoid wrong resetting when not yet returned
            end
            % STE: was else if
            if (agents(m,n) >= ideas_space_size)
                %agents(m,n)=ideas_space_size;

                % Important Difference!!! whether Bouncing or Zero v
                %v(m,n)=0;
                v(m,n)=-abs(v(m,n));
            end
        end
    end
   
     
    %% Euler Step
    tattr = ths(agents);
    leap = v + force + tattr;
    
    if (params.noisetype == noise_on_p)
        %Add noise - it was already added if it was on velocity
        leap = leap + normrnd(0,sigma,ideas_space_dim,n_agents);
    end
    
    agents = agents + leap.*dt;
    
   
    %agents;
    
    % NaN values can be produced due to limited machine precision
    % if agents are very close to pos. 0
    % in this case, we can just reset them to 0
    for i=1:ideas_space_dim
        problematic = isnan(agents(i,:));
        agents(i,problematic) = 0;
    end
    
    %% Statistics, Plot and Dump
    if (VIDEO)
        
        switch (params.plottype)
        
            case plot_cross
            % PLOT red crosses
            plot(agents(1,:),agents(2,:),'rx');     
        
            case plot_number
            % PLOT BLACK NUMBERS
            text(agents(1,:),agents(2,:), num2str([1:length(agents)]'));
            plot(exp(agents(1,1)));
        
            case plot_number_color
            % PLOT COLORED NUMBERS
            points = arrayfun(@(x) {[ '\color{' colors{mod(x,length(colors))+1} '}' int2str(x)]}, 1:length(agents));
            text(agents(1,:),agents(2,:), points');

        end
        hold on;
        plot(truth(1),truth(2),'go');
        hold off;
            
        xlim([0 ideas_space_size]);
        ylim([0 ideas_space_size]);
        
            
        if (DEBUG)
            %Debugging: Calculate energy and control whether it is constant
            E=0;
            for i=1:n_agents
                E=E+0.5*(v(:,i))'*v(:,i);
            end
            energy=num2str(E);
            legend(energy);
        end
        
        pause(0.1);
    end
    
    
    if (DUMP)
        if ( mod(t / dt, DUMP_RATE) == 0)
            agentspos = cat(3,agentspos,agents);
            agentsv = cat(3,agentsv,v);
        end
    end
      
        
   
end         %end of time loop

%% Calculation of a measure for the convergence to the TRUTH

%conv=norm(agents-repmat(truth,1,n_agents));    %possibly wrong
elements = zeros(n_agents,1);
for j=1:n_agents
   elements(j) = norm(agents(:,j)-truth); 
end

conv=sqrt((elements'*elements)/n_agents); %Root Mean Square (=Quadratic Mean)

%% Dump simulation data
if (DUMP)
    
    % Last position is added anyway
    if ( mod(t / dt, DUMP_RATE) ~= 0)
        agentspos = cat(3,agentspos,agents);
        agentsv = cat(3,agentsv,v);
     end
    
    
    
    
    
    %dump = {parameters,agentsdump,truth,conv};
    dump = struct('name', folderName, ...
                  'sim', simCount, ...
                  'run', nRun, ...
                  'parameters', params,...
                  'agents', agentspos, ...
                  'agentsv', agentsv, ...
                  'truth',truth,...
                  'conv', conv...
                 );
    
                
    [status,message,messageid] = mkdir(folderName);               
                
    fileName = [folderName '/' int2str(simCount) '-' int2str(nRun) '.mat']; 
    
    save(fileName,'dump');
    
end
 
toc(tStart);
end
