function param_sets_LSF (params)

parallel.importProfile('/cluster/apps/matlab/support/BrutusLSF8h.settings')

TASKS4JOB = 20; % How many tasks group in one job
jobCount = 1;

logFolder = ['log/' params.simName];
mkdir(logFolder); % The name is unique under the dump directory.
dumpFolder = [ params.dumpDir params.simName];

sched = findResource('scheduler','type','lsf');
sched = findResource('scheduler','type','BrutusLSF8h');
emails = 'sbalietti@ethz.ch';
%submitArgs = ['-o ' logFolder '/' simName '.log -u ' emails];
%submitArgs = ['-o ' logFolder '/' simName '.log -B']; -B / -N sends email
submitArgs = ['-o ' logFolder '/' params.simName '.log'];
set(sched, 'SubmitArguments',submitArgs);
set(sched, 'DataLocation', [logFolder '/']);


jobName = genvarname(['j_' int2str(jobCount)]);
eval([jobName '= createJob(sched)']);
jobCount = jobCount + 1;

nCombinations = size(params.dts,2)*size(params.n_agents,2)*size(params.ideas_space_sizes,2)*...
                size(params.ideas_space_dims,2)*size(params.As,2)*size(params.Bs,2)*size(params.ks,2)*...
                size(params.d0s,2)*size(params.d1s,2)*size(params.alphas,2)*size(params.taus,2)*size(params.Rs,2)*...
                size(params.sigmas,2)*size(params.v_scalings,2)*size(params.nof_clusters,2)*...
                size(params.clusterTightness,2)*size(params.truths,2);
      
simCount=1; %counter of all simulations
%nest several loops to simulate parameter sets
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
        
        for rCount=1:params.nRuns
        
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
                'attrtype', attrtype ...
            );
            
            createTask(eval(jobName), @simulation, 0, {paramsObj});


           % Submit the job to the scheduler if a sufficient number of task
           % is reached
           if (mod(simCount,TASKS4JOB)==0)
               submit(eval(jobName))

               if (simCount ~= nCombinations)
                   jobName = genvarname(['j_' int2str(jobCount)]);
                   eval([jobName '= createJob(sched)']);
                   jobCount = jobCount + 1;
                   %[jobName,jobCount] = createJobName(jobCount);
               end

           end

            fprintf('\n\n');
        end
        simCount=simCount+1; %updating the simulations count
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

% Submit the left-over tasks
if (mod(simCount,TASKS4JOB) ~= 1) % 1: it has always one last increment
    submit(eval(jobName))
end


end

% function [jobName,jobCount] = createJobName(jobCount) 
%    jobName = genvarname(['j_' int2str(jobCount)]);
%    eval([jobName '= createJob(sched)']);
%    jobCount = jobCount + 1;
% end



     