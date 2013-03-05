function param_sets_local(dumpDir,simName,VIDEO,DEBUG,DUMP,...
    nRuns, dts,t_ends,n_agents,ideas_space_sizes,ideas_space_dims,...
    As,Bs,ks,d0s,d1s,alphas,taus,Rs,sigmas,...
    v_scalings,nof_clusters,cts,truths)


dumpFolder = [dumpDir simName];


nCombinations = size(dts,2)*size(n_agents,2)*size(ideas_space_sizes,2)*...
                size(ideas_space_dims,2)*size(As,2)*size(Bs,2)*size(ks,2)*...
                size(d0s,2)*size(d1s,2)*size(alphas,2)*size(taus,2)*size(Rs,2)*...
                size(sigmas,2)*size(v_scalings,2)*size(nof_clusters,2)*...
                size(cts,2)*size(truths,2);

%nest several loops to simulate parameter sets


simCount=1; % counter of all simulations within a Run

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
        truth = truths(:,i16);
        
        for rCount=1:nRuns
        

            fprintf('\n%s\n',simName);
            fprintf('Starting Run: %d/%d of Simulation n. %d/%d:\n', ...
                    rCount,nRuns,simCount,nCombinations)
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
            fprintf('%+15s = %2.3f\n','v scaling',v_scaling);
            fprintf('%+15s = %1d\n','n cluster',nof_cluster);
            fprintf('%+15s = %2.3f\n','tight clusters',clusterTightness);
            fprintf('%+15s = [%2.3f:%2.3f]\n','truth',truth(1,1),truth(2,1));
            fprintf('------------------------------------\n');

            simulation(dumpFolder,simCount,VIDEO,DEBUG,DUMP,...
                                      rCount,dt,t_end,agents,iss,isd,...
                                      A,B,k,d0,d1,alpha,tau,R,sigma,...
                                      v_scaling,nof_cluster,...
                                      clusterTightness,truth);

            fprintf('\n\n');
            
        end % End n runs with identical param set
             
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