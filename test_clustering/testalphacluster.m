%% Add other directories to path
path(path,'../util/'); % Help functions
path(path,'../lib/'); % Help functions
path(path,'../'); % Help functions

%params

VIDEO = 1;
SHOW_POTENTIAL = 0;

DUMP = 0;
DUMP_RATE = 1;
nRun = 1;
dt = 0.01;
t_end = 2;
n_agents = 10;
ideas_space_size = 1;
ideas_space_dim = 2;

A = 2;
B = 2;
d0 = 1;
d1 = 1;

alphaParam = 0.1;
tau = 0.1;
R = 0.03;
sigma = 0.01;
epsilon = 0.1;

truth = [0.4; 0.4];

LIMS = [0.3, 0.6];

noisetype = 4;
boundaries = 1;

SEED = 0;

display('Starting...');


s = RandStream('mcg16807','Seed', SEED);
    RandStream.setGlobalStream(s);
    
%% Initializations and Definitions


% Row-vector holding agents' positions.
agents = 0.5 .* ones(2, n_agents);

v = 0.5 .* ones(2, n_agents-1);
v(:,n_agents) = [0;0.5];

rownorm = @(X,P) sum(abs(X).^P,2).^(1/P);
colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);
    
truth4all = repmat(truth, 1, n_agents);        
    
plottype = 3;

% PLOT TYPE
plot_cross = 0;
plot_number = 1;
plot_number_color = 2;
plot_arrow = 3;


% Attraction to truth force type
attr_zero = 0;
attr_const = 1;
attr_linear = 2;
attr_expo = 3;
attr_millean_arena = 4;
attr_hard_to_find = 5;
attr_wide_funnel = 6;
attr_gentle_landing = 7;

ths = @(x) (truth4all - x) ./ tau;
    

% NOISE TYPES
noise_on_p = 0;
noise_on_v = 1;
noise_adaptive_on_v =2;
noise_on_v_angular = 3;
noise_navp = 4;



% d: Eucledian distance points i,j
d = @(agents,i,j) norm(agents(:,i)-agents(:,j));  

% df: Attraction force between agents i and j
df = @(agents,i,j) ( A*d(agents,i,j)^k * exp( -d(agents,i,j) / d0) - B*d(agents,i,j)^k * exp( -d(agents,i,j) / d1)) .* ( (agents(:,j) - agents(:,i)) / d(agents,j,i) );                 

%% Initial Plot
if (VIDEO)
    fig1 = figure(1);

    hold on
    
    % Showing the potential of the attraction to truth
    if (SHOW_POTENTIAL)
        PRECISION = 0.01;
        a = [0:PRECISION:ideas_space_size;0:PRECISION:ideas_space_size];
        [X,Y] = meshgrid(a(1,:), a(1,:));
        Z = zeros(size(X));
        for i=1:length(X)
            potential_grid = [X(i,:) ; Y(i,:)];
            forces = ths(potential_grid);
            Z(i,:) = colnorm(forces,2);
        end
    end
    
    
    switch (plottype)

        case plot_cross
        % PLOT red crosses
        % plot(agents(1,:),agents(2,:),'rx');     
        plot(agents(1,:), agents(2,:),'ok', ...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','r',...
            'MarkerSize', 5);

        case plot_number
        % PLOT BLACK NUMBERS
        text(agents(1,:),agents(2,:), num2str([1:length(agents)]'));
        plot(exp(agents(1,1)));

        case plot_number_color
        % PLOT COLORED NUMBERS
        cla;
        points = arrayfun(@(x) {[ '\color{' colors{mod(x,length(colors))+1} '}' int2str(x)]}, 1:length(agents));
        text(agents(1,:),agents(2,:), points');

        case plot_arrow
        % PLOT VELOCITY ARROWS
        quiver(agents(1,:),agents(2,:),v(1,:),v(2,:));


    end
    hold on;
        
    if (SHOW_POTENTIAL)
        [C, h] = contour(X,Y,Z);
        alpha(.5);
    end
   
    
    plot(truth(1),truth(2),'go');
    
    hold off
    
    xlim(LIMS);
    ylim(LIMS);
end

            
        

%if (DUMP)
    agentspos = agents;
    agentsv = v;
%end

%% Time iteration

% Start Timer
tStart = tic; 

counter = 0;
for t=0:dt:t_end 
  
noConsensus = 1;

% while( noConsensus)

    counter = counter + 1;
    
    %% Compute Total Attraction Force
    
    % agents are columns (i), dims stored in rows, sum every entry over j
    force = zeros(ideas_space_dim,n_agents);                     
    
    %% Compute Velocity influenced by agents in cutoff
    
    
    avgvs = zeros(n_agents,1);
    
    for i=1:n_agents
        % Important: the update is sequential and random in the order
        other_agents_ids = mod([i:n_agents+i-2],n_agents)+1;
        other_agents_ids = other_agents_ids(randperm(length(other_agents_ids)));
        agents_in_cutoff = zeros(1,n_agents);
        for ii=1:size(other_agents_ids,2)
            j = other_agents_ids(ii);
            if (d(agents,i,j) <= R)
                agents_in_cutoff(j) = 1;
                %v(:,i)=v(:,i)+(1-alpha).*av_velocity(i,R); 
                %agents_in_cutoff=find(d(i,j)<R);
            end
        end

        % Important: only if some agents are in the cutoff radius
        % the velocity is split between alpha and 1-alpha.
        % If not the agent continuous his search undisturbed
        if (nnz(agents_in_cutoff) ~=0) 
            

                average_velocity=(v*(agents_in_cutoff)') ./ nnz(agents_in_cutoff);
                v(:,i) = alphaParam * v(:,i) + (1-alphaParam) * average_velocity;
                avgvs(i) = norm(average_velocity);

        end
       
        if (noisetype == noise_on_v_angular || noisetype == noise_navp)
            % Angular noise
            % we lose something in the conversion. But it is not increasing
            % over time (maybe -- last check was OK).
            [THETA, RHO] = cart2pol(v(1,i), v(2,i));
            THETA2 = normrnd(THETA, sigma);
            [vX, vY] = pol2cart(THETA2, RHO);            
            v(1,i) = vX;
            v(2,i) = vY;            
            
        elseif (noisetype == noise_on_v)
            % Important: noise on velocity
            %v(:,i) = v(:,i) + normrnd(0, sigma, ideas_space_dim, 1);
            v(:,i) = normrnd(v(:,i), sigma, ideas_space_dim, 1);
            
        elseif (noisetype == noise_adaptive_on_v)
            % Adaptive noise: the faster you go the higher the noise
            v(:,i) = normrnd(v(:,i), sigma, ideas_space_dim, 1);    
        end
            
        
    end    
    
    %% Computing attraction force
    tattr = ths(agents);
    
    if (boundaries == 1)
    %% Boundary Conditions: reflecting from the walls
        for n = 1 : n_agents        
           for m = 1 : ideas_space_dim 

                if (agents(m,n) <= 0)
                    v(m,n) = abs(v(m,n));
                    % v(m,n) = abs(v(m,n)) - tattr(m,n);
                elseif (agents(m,n) >= ideas_space_size)
                    v(m,n) = - abs(v(m,n));
                    %v(m,n) = - abs(v(m,n)) + tattr(m,n);

                end
          end

        end
    end   
     
    %% Euler Step

   
    % Forces and v stay separated
    leap = v + tattr + force;

    
    if (noisetype == noise_navp || noisetype == noise_on_p)
        % Add noise
        leap = leap + normrnd(0, epsilon, ideas_space_dim, n_agents);
    end
    
    % mean(colnorm(v,2))
    
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
        
        switch (plottype)
        
            case plot_cross
            % PLOT red crosses
            % plot(agents(1,:),agents(2,:),'rx');     
            plot(agents(1,:), agents(2,:),'ok', ...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize', 5);
            
            case plot_number
            % PLOT BLACK NUMBERS
            text(agents(1,:),agents(2,:), num2str([1:length(agents)]'));
            plot(exp(agents(1,1)));
        
            case plot_number_color
            % PLOT COLORED NUMBERS
            cla;
            points = arrayfun(@(x) {[ '\color{' colors{mod(x,length(colors))+1} '}' int2str(x)]}, 1:length(agents));
            text(agents(1,:),agents(2,:), points');
           
            case plot_arrow
            % PLOT VELOCITY ARROWS
            quiver(agents(1,:),agents(2,:),v(1,:),v(2,:));
            

        end
        hold on;
        
        if (SHOW_POTENTIAL)
            [C, h] = contour(X,Y,Z);
            alpha(.5);
        end
        
        plot(truth(1),truth(2),'go');
        
        hold off;
            
        if (boundaries ~= 0)
            xlim(LIMS);
            ylim(LIMS);
        end
        
        pause(1);
        
    end
    
    if (DUMP)
        if ( mod(counter, DUMP_RATE) == 0)
            agentspos = cat(3,agentspos,agents);
            agentsv = cat(3,agentsv,v);
        else
            %t
            %counter
            %mod(t / dt, DUMP_RATE)
        end
    end
    
    % Stop condition
%     norm_from_truth = colnorm(agents - truth4all, 2);
%     mean_agents_fromtruth = mean(norm_from_truth);
% 
%     if (mean_agents_fromtruth < 0.05)
%         noConsensus = 0;
%     elseif (counter > 20000)
%         counter
%         noConsensus = 0;
%     end
    
end
% end of time loop.


%% Dump simulation data
if (DUMP)
    
    dump = struct('name', folderName, ...
                  'sim', simCount, ...
                  'run', nRun, ...
                  'agents', agentspos, ...
                  'agentsv', agentsv, ...
                  'truth',truth, ...                 
                  'counter', counter ...
                 );
    
                
    [status,message,messageid] = mkdir(folderName);               
                
    fileName = [folderName '/' int2str(simCount) '-' int2str(nRun) '.mat']; 
    
    save(fileName,'dump');
    
end
toc(tStart);
