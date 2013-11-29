function param_sets_LSF (params)

% SEED TYPE
seed_fixed = 0;
seed_random = 1;
seed_machinetime = 2;

if (params.seedtype ~= seed_fixed)
    s = RandStream('mcg16807','Seed', params.batchSeed);
    RandStream.setGlobalStream(s);
end 

path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

parallel.importProfile('/cluster/apps/matlab/support/BrutusLSF8h.settings')


% How many sequential simulations in one task.
SIMS4TASK = 4; 
% How many tasks group in one job. Should be multiple with SIMS4TASK.
TASKS4JOB = 20;

jobCount = 1;

logFolder = ['log/' params.simName];
mkdir(logFolder); % The name is unique under the dump directory.
dumpFolder = [ params.dumpDir params.simName];

sched = findResource('scheduler','type','lsf');
%sched=parcluster('BrutusLSF8h');
emails = 'sbalietti@ethz.ch';
%submitArgs = ['-o ' logFolder '/' simName '.log -u ' emails];
%submitArgs = ['-o ' logFolder '/' simName '.log -B']; -B / -N sends email
submitArgs = [' -R "rusage[mem=8000]" -o ' logFolder '/' params.simName '.log'];
set(sched, 'SubmitArguments',submitArgs);
set(sched, 'DataLocation', [logFolder '/']);

j = createJob(sched);

jobCount = jobCount + 1;

nCombinations = size(params.dts,2)*size(params.n_agents,2)*size(params.ideas_space_sizes,2)*...
                size(params.ideas_space_dims,2)*size(params.As,2)*size(params.Bs,2)*size(params.ks,2)*...
                size(params.d0s,2)*size(params.d1s,2)*size(params.alphas,2)*size(params.taus,2)*size(params.Rs,2)*...
                size(params.sigmas,2)*size(params.v_scalings,2)*size(params.nof_clusters,2)*...
                size(params.clusterTightness,2)*size(params.truths,2)*size(params.forces_on_v,2);

% Init cell array of cell arrays           
paramObjs = cell(SIMS4TASK,1);            

% Counter of all simulations.
simCount = 1;

% Nest several loops to simulate parameter sets.
for i1=1:size(params.dts)
    dt = params.dts(i1);
     
    for i2=1:size(params.t_ends,2)
        t_end = params.t_ends(i2);   
    
    for i3=1:size(params.n_agents,2)
        agents = params.n_agents(i3);
        
    for i4=1:size(params.ideas_space_sizes,2)
        iss = params.ideas_space_sizes(i4);
                
    for i5=1:size(params.ideas_space_dims,2)
        isd = params.ideas_space_dims(i5);
         
    for i6=1:size(params.As,2)
        A = params.As(i6);     
        
    for i6b=1:size(params.Bs,2)
        B = params.Bs(i6b);
        
    for i7=1:size(params.ks,2)
        k = params.ks(i7);
        
    for i8=1:size(params.d0s,2)
        d0 = params.d0s(i8);
        
    for i8b=1:size(params.d1s,2)
        d1 = params.d1s(i8b);
        
    for i9=1:size(params.alphas,2)
        alpha = params.alphas(i9);
        
    for i10=1:size(params.taus,2)
        tau = params.taus(i10);    
           
    for i11=1:size(params.Rs,2)
        R = params.Rs(i11);    
        
    for i12=1:size(params.sigmas,2)
        sigma = params.sigmas(i12); 
        
    for i13=1:size(params.v_scalings,2)
        v_scaling = params.v_scalings(i13);
    
    for i14=1:size(params.nof_clusters,2)
        nof_cluster = params.nof_clusters(i14);
        
    for i15=1:size(params.clusterTightness,2)
        clusterTightness = params.clusterTightness(i15);    
        
    for i16=1:size(params.truths,2)
        truth = params.truths(:,i16);
        
    for i17=1:size(params.attrtype,2)
        attrtype = params.attrtype(:,i17);
        
    for i18=1:size(params.noisetype,2)
        noisetype = params.noisetype(:,i18);
        
    for i19=1:size(params.forces_on_v,2)
        forces_on_v = params.forces_on_v(:,i19);
        
        for rCount=1:params.nRuns
            
            % Defining seed
            if (params.seedtype == seed_machinetime)
                % Set seed with milliseconds
                seed1 = sscanf(datestr(now, 'FFF'),'%d') * 1000;
                s = RandStream('mcg16807','Seed', seed1);
                RandStream.setGlobalStream(s);
                rng shuffle
                seed2 = randi(1000);
                seed = seed1 + seed2;
                % Random pause to break correlation in itialization of
                % simulations
                pause(randi(1));
            elseif (params.seedtype == seed_random)
                % Random seed
                seed = randi(1000000);
                
            elseif (params.seedtype == seed_fixed)
                seed = params.seed;
            end
        
            fprintf('\n%s\n',params.simName);
            fprintf('Starting Run: %d/%d of Simulation n. %d/%d:\n', ...
                     rCount,params.nRuns,simCount,nCombinations)
            fprintf('------------------------------------\n');
            fprintf('%+15s = %d\n','steps',t_end);
            fprintf('%+15s = %2f\n','dt',dt);
            fprintf('%+15s = %d\n','nAgents',agents);  
            fprintf('%+15s = %1d\n','IdeasSpace size',iss);
            fprintf('%+15s = %1d\n','IdeasSpace dims',isd);    
            fprintf('%+15s = %2.3f\n','A',A); 
            fprintf('%+15s = %2.3f\n','B',B); 
            fprintf('%+15s = %2.3f\n','k',k);  
            fprintf('%+15s = %2.3f\n','d0',d0);
            fprintf('%+15s = %2.3f\n','d1',d1);
            fprintf('%+15s = %2.3f\n','alpha',alpha);
            fprintf('%+15s = %2.3f\n','tau',tau);   
            fprintf('%+15s = %2.3f\n','R',R);
            fprintf('%+15s = %2.3f\n','sigma',sigma);
            fprintf('%+15s = %2.3f\n','v_scaling',v_scaling);
            fprintf('%+15s = %1d\n','n cluster',nof_cluster);
            fprintf('%+15s = %2.3f\n','tight clusters',clusterTightness);
            fprintf('%+15s = [%2.3f:%2.3f]\n','truth',truth(1,1),truth(2,1));
            fprintf('%+15s = %d\n', 'Attr. type', attrtype);
            fprintf('%+15s = %d\n', 'Noise type', noisetype);
            fprintf('%+15s = %d\n', 'Plot type', params.plottype);
            fprintf('%+15s = %d\n', 'Forces on V', forces_on_v);
            fprintf('%+15s = %d\n','Seed', seed);
            fprintf('------------------------------------\n');

            paramsObj = struct( ...
                'folderName', dumpFolder, ...
                'simCount', simCount, ...
                'VIDEO', params.VIDEO, ...
                'DEBUG', params.DEBUG, ...
                'DUMP', params.DUMP, ...
                'DUMP_RATE', params.DUMP_RATE, ...
                'seedtype', params.seedtype, ...
                'plottype', params.plottype, ...
                'SHOW_POTENTIAL', params.SHOW_POTENTIAL, ...
                'nRun', rCount, ...
                'dt', dt, ...
                't_end', t_end, ...
                'n_agents', agents, ...
                'ideas_space_size', iss, ...
                'ideas_space_dim', isd, ...
                'A', A, ...
                'B', B, ...
                'k', k, ...
                'd0', d0, ...
                'd1', d1, ...
                'alpha', alpha, ...
                'tau', tau, ...
                'R', R, ...
                'sigma', sigma, ...
                'v_scaling', v_scaling, ...
                'nof_cluster', nof_cluster, ...
                'clusterTightness', clusterTightness, ...
                'truth', truth, ...
                'noisetype', noisetype, ...
                'attrtype', attrtype, ...
                'forces_on_v', forces_on_v, ...
                'seed', seed ...
            );
  
            % It is convenient to group together more simulations in one
            % task if simulations are short. Matlab overhead to start on
            % each cluster node is about 1 minute.
            taskIdx = mod(simCount, SIMS4TASK);

            if (taskIdx ~= 0)
                paramObjs{taskIdx} = {paramsObj};
            else
                paramObjs{taskIdx+1} = {paramsObj};
                createTask(j, @simulation, 0, paramObjs);
                paramObjs = cell(SIMS4TASK, 1);
            end
  
            % createTask(j, @simulation, 0, {paramsObj});

            % Submit the job to the scheduler in batches
            if (mod(simCount, TASKS4JOB)==0)
                submit(j);

                if (simCount ~= nCombinations)
                    j = createJob(sched); 
                    jobCount = jobCount + 1;
                end

            end

            fprintf('\n\n');
        end
        % updating the simulations count
        simCount = simCount + 1;
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
end  

% Submit the left-over tasks.
if (mod(simCount, TASKS4JOB) ~= 1) % 1: it has always one last increment.
    if (taskIdx ~= 0)
        createTask(j, @simulation, 0, paramObjs);
    end
    submit(j);
end


end



     
