function param_sets_parallel(dumpDir,simName,VIDEO,DEBUG,DUMP,...
    nRuns,dts,t_ends,n_agents,ideas_space_sizes,ideas_space_dims,...
    As,Bs,ks,d0s,d1s,alphas,taus,Rs,sigmas,...
    v_scalings,nof_clusters,cts,truths)


% IMPORTANT: Not Working Weird Error on ParFor

dumpFolder = [ dumpDir simName];

simCount=1; %counter of all simulations

nCombinations = size(dts,2)*size(n_agents,2)*size(ideas_space_sizes,2)*...
                size(ideas_space_dims,2)*size(As,2)*size(Bs,2)*size(ks,2)*...
                size(d0s,2)*size(d1s,2)*size(alphas,2)*size(taus,2)*size(Rs,2)*...
                size(sigmas,2)*size(v_scalings,2)*size(nof_clusters,2)*...
                size(cts,2)*size(truths,2);
            
% Notice: truths and ideas_space_dims must be consistent over all the param
% range.
            
% Init Sliced Arrays of params
dtss = zeros(1,nCombinations);
t_endss = zeros(1,nCombinations);
agentss = zeros(1,nCombinations);
issss = zeros(1,nCombinations);
isdss = zeros(1,nCombinations);
Ass = zeros(1,nCombinations);
Bss = zeros(1,nCombinations);
kss = zeros(1,nCombinations);
d0ss = zeros(1,nCombinations);
d1ss = zeros(1,nCombinations);
alphass = zeros(1,nCombinations);
tauss = zeros(1,nCombinations);
Rss = zeros(1,nCombinations);
sigmass = zeros(1,nCombinations);
v_scalingss = zeros(1,nCombinations);
nof_clusterss = zeros(1,nCombinations);
ctss = zeros(1,nCombinations);
truthss = zeros(size(truths,1),nCombinations);

% TODO: there is probably a method to do this, like when you do the
% meshgrids. Investigate it.

%nest several loops to simulate parameter sets
for i1=1:size(dts,2)
    dt = dts(i1);
       
    for i2=1:size(t_ends,2)
        t_end = t_ends(i2);
          
    for i3=1:size(n_agents,2)
        agents = n_agents(i3);
                 
    for i4=1:size(ideas_space_sizes,2)
        iss = ideas_space_sizes(i4);
                
    for i5=1:size(ideas_space_dims,2)
        isd = ideas_space_dims(i5);
                 
    for i6=1:size(As,2)
        A = As(i6);     
                
    for i6b=1:size(Bs,2)
        B = Bs(i6b);
                
    for i7=1:size(ks,2)
        k = ks(i7);
                   
    for i8=1:size(d0s,2)
        d0 = d0s(i8);
                
    for i8b=1:size(d1s,2)
        d1 = d1s(i8b);
            
    for i9=1:size(alphas,2)
        alpha = alphas(i9);
            
    for i10=1:size(taus,2)
        tau = taus(i10);    
             
    for i11=1:size(Rs,2)
        R = Rs(i11);    
            
    for i12=1:size(sigmas,2)
        sigma = sigmas(i12); 
        
    for i13=1:size(v_scalings,2)
        v_scaling = v_scalings(i13);
    
    for i14=1:size(nof_clusters,2)
        nof_cluster = nof_clusters(i14);
        
    for i15=1:size(cts,2)
        clusterTightness = cts(i15);
        
    for i16=1:size(truths,2)
        truth = truths(i16);
        
                
        dtss(simCount) = dt;
        t_endss(simCount) = t_end; 
        agentss(simCount) = agents;
        issss(simCount) = iss;
        isdss(simCount) = isd;
        Ass(simCount) = A;
        Bss(simCount) = B;
        kss(simCount) = k;
        d0ss(simCount) = d0;
        d1ss(simCount) = d1;
        alphass(simCount) = alpha;
        tauss(simCount) = tau;
        Rss(simCount) = R;
        sigmass(simCount) = sigma;
        v_scalingss(simCount) = v_scaling;
        nof_clusterss(simCount) = nof_cluster;
        ctss(simCount) = clusterTightness;       
        truthss(:,simCount) = truth;

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

matlabpool open

% Notice parfor cannot be nested. You should parallelize the loop with 
% more iteration between the Simulation Loop and the Runs Loop

parfor i=1:nCombinations % Simulations Loop
    
     dt = dtss(i);
     t_end = t_ends(i);
     agents = agentss(i);
     iss = issss(i);
     isd = isdss(i);
     A = Ass(i);
     B = Bss(i);
     k = kss(i);
     d0 = d0ss(i);
     d1 = d1ss(i);
     alpha = alphass(i);
     tau = tauss(i);
     R = Rss(i);
     sigma = sigmass(i);
     v_scaling = v_scalingss(i);
     nof_cluster = nof_clusterss(i);
     clusterTightness = ctss(i);
     truth = truthss(:,i);
    
     
    for rCount=1:nRuns % Runs Loop
        

        fprintf('\n%s\n',simName);
        fprintf('Starting Run: %d/%d of Simulation n. %d/%d:\n', ...
                rCount,nRuns,i,nCombinations)
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
        fprintf('------------------------------------\n');


        simulation(dumpFolder,i,VIDEO,DEBUG,DUMP,...
                    rCount,dt,t_end,agents,iss,isd,...
                    A,B,k,d0,d1,alpha,tau,R,sigma,...
                    v_scaling,nof_cluster,...
                    clusterTightness,truth);

    end
    
end

matlabpool close
% Notice: If an error occurred the pool will not be closed, and must be closed
% manually