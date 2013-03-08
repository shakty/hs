function [row] = create_params_string(params, truth) 

    row = [ ' dt: ' num2str(params.dt) ', steps: ' num2str(params.steps) ', agents: ' num2str(params.nAgents)];
    
    row = [ row ', size: ' num2str(params.iss) ',alpha: ' num2str(params.alpha)];
    
    row = [ row ', R: ' num2str(params.R) ', k: ' num2str(params.k) ', A:' num2str(params.A)];
    
    row = [ row ', d0: ' num2str(params.d0) ', B: ' num2str(params.B) ', d1: ' num2str(params.d1)];
    
    row = [ row ', tau: ' num2str(params.tau) ', sigma: ' num2str(params.sigma) ', v0: ' num2str(params.vScaling)];
    
    row = [ row ', c0: ' num2str(params.nClusters) ',ct0: ' num2str(params.clusterTightness)];
    
    row = [ row ', truth: ' num2str(truth(1,1)) ';' num2str(truth(2,1))];
end