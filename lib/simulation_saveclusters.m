function conv = simulation_saveclusters (folderName,simCount,VIDEO,DEBUG,DUMP,...
    nRun, dt,t_end,n_agents,ideas_space_size,ideas_space_dim,...
    A,B,k,d0,d1,alpha,tau,R,sigma,v_scaling,nof_cluster,clusterTightness,...
    truth)


display('Starting...');

%% Initializations and Definitions

% Row-vector holding agents' positions.
% agents = ideas_space_size.*rand(ideas_space_dim,n_agents);
agents = initial_pos_clustered(nof_cluster,clusterTightness,n_agents,...
                                ideas_space_size,ideas_space_dim);

% center of ideas
agents_average=zeros(ideas_space_dim,1); 


% Row-vector holding agents' velocities
v = v_scaling*(rand(ideas_space_dim,n_agents)-0.5);

%function handle to ODE-rhs: v + attraction force + noise + convergenge to truth
rhs = @(x,v,force) v + force + normrnd(0,sigma,ideas_space_dim,n_agents) + (repmat(truth,1,n_agents)-x)./tau;                       

% d: Eucledian distance points i,j
d = @(agents,i,j) norm(agents(:,i)-agents(:,j));  

% df: Attraction force between agents i and j
df = @(agents,i,j) ( A*d(agents,i,j)^k * exp( -d(agents,i,j) / d0) - B*d(agents,i,j)^k * exp( -d(agents,i,j) / d1)) .* ( (agents(:,j) - agents(:,i)) / d(agents,j,i) );                 

%% Initial Plot
if (VIDEO)
    fig1 = figure(1);

    % Initial Positioning
    %hold on
    plot(agents(1,:),agents(2,:),'rx');
    %plot(agents_average(1),agents_average(2),'bx');
    plot(truth(1),truth(2),'go');
end


if (DUMP)
    agentsdump = agents;
end

%% Time iteration

% Start Timer
tStart = tic; 


for t=0:dt:t_end                      
    % update state (Euler)
    
    %% Average Position
    for i=1:ideas_space_dim
        agents_average(i)=sum(agents(i,:))/n_agents;
    end
    
    %% Compute Total Attraction Force
    
    % agents are columns (i), dims stored in rows, sum every entry over j
    force = zeros(ideas_space_dim,n_agents);                     
    for i=1:n_agents
        % Generates the ids of all the other agents
        other_agents_ids = mod([i:n_agents+i-2],n_agents)+1;
        for ii=1:size(other_agents_ids)
            j = other_agents_ids(ii);
            force_current = df(agents,i,j);
            force(:,i)= force(:,i) + force_current;
        end
     end
    
    %% Compute Velocity influenced by agents in cutoff  
    for i=1:n_agents
        other_agents_ids = mod([i:n_agents+i-2],n_agents)+1;
        agents_in_cutoff=zeros(1,n_agents);
        for ii=1:size(other_agents_ids,2)
            j = other_agents_ids(ii);
            if (d(agents,i,j) <= R)
                agents_in_cutoff(j)=1;
                %v(:,i)=v(:,i)+(1-alpha).*av_velocity(i,R); 
                %agents_in_cutoff=find(d(i,j)<R);
            end
        end

        if(nnz(agents_in_cutoff)==0) %interpret average velocity as 0 if noone in cutoff
           average_velocity=0;
        else
            average_velocity=(v*(agents_in_cutoff)')./nnz(agents_in_cutoff);
        end
        
        v(:,i)= alpha*v(:,i)+(1-alpha)*average_velocity;  
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
    agents = agents + rhs(agents,v,force).*dt;
    
    % NaN values can be produced due to limited machine precision
    % if agents are very close to pos. 0
    % in this case, we can just reset them to 0
    for i=1:ideas_space_dim
        problematic = isnan(agents(i,:));
        agents(i,problematic) = 0;
    end
    
    %% Statistics, Plot and Dump
    if (VIDEO)
        
        %plot(exp(agents(1,1)));
        plot(agents(1,:),agents(2,:),'rx'); 
        hold on;
        plot(agents_average(1),agents_average(2),'bo');
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
        
        pause(0.01);
    end
    
    
    if (DUMP)
        agentsdump = cat(3,agentsdump,agents);
    end
      
    %print(gcf, '-dpng', 'testPlot.png', '-r 10')
        
   
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
    
    
    %parameters = [n_agents,A,B,k,d0,d1,alpha,tau,R,sigma,v_scaling,nof_cluster];
    parameters = struct(...
                        'steps', t_end, ...
                        'dt', dt, ...
                        'nAgents', n_agents, ...
                        'iss', ideas_space_size, ...
                        'isd', ideas_space_dim, ...
                        'A', A, ...
                        'B', B, ...
                        'k',k, ...
                        'd0', d0, ...
                        'd1',d1, ...
                        'alpha',alpha, ...
                        'tau',tau, ...
                        'R',R, ...
                        'sigma',sigma, ...
                        'vScaling', v_scaling, ...
                        'nClusters',nof_cluster, ...
                        'clusterTightness',clusterTightness ...
                        );     
    
    
    %dump = {parameters,agentsdump,truth,conv};
    dump = struct('name', folderName, ...
                  'sim', simCount, ...
                  'run', nRun, ...
                  'parameters',parameters,...
                  'agents',agentsdump, ...
                  'truth',truth,...
                  'conv', conv...
                 );
    
                
    [status,message,messageid] = mkdir(folderName);               
                
    fileName = [folderName '/' int2str(simCount) '-' int2str(nRun) '.mat']; 
    
    save(fileName,'dump');
    
end

 
toc(tStart);
end