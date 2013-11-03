function [row] = create_params_string(params, truth) 

    row = [ ' dt: ' num2str(params.dt) ', steps: ' num2str(params.t_end) ', agents: ' num2str(params.n_agents)];
    
    row = [ row ', size: ' num2str(params.ideas_space_size) ',alpha: ' num2str(params.alpha)];
    
    row = [ row ', R: ' num2str(params.R) ', k: ' num2str(params.k) ', A:' num2str(params.A)];
    
    row = [ row ', d0: ' num2str(params.d0) ', B: ' num2str(params.B) ', d1: ' num2str(params.d1)];
    
    row = [ row ', tau: ' num2str(params.tau) ', sigma: ' num2str(params.sigma) ', v0: ' num2str(params.v_scaling)];
    
    row = [ row ', c0: ' num2str(params.nof_cluster) ',ct0: ' num2str(params.clusterTightness)];
    
    row = [ row ', truth: ' num2str(truth(1,1)) ';' num2str(truth(2,1))];
end