function [row] = create_params_string_small(params) 

    row = [ 'alpha: ' num2str(params.alpha) ' R: ' num2str(params.R)];
    
    row = [ row ' sigma: ' num2str(params.sigma)];

    row = [ row ' agents: ' num2str(params.n_agents)];
        
    row = [ row ' tau: ' num2str(params.tau)];
    
    row = [ row ' v0: ' num2str(params.v_scaling)];
    
    row = [ row ' noise: ' num2str(params.noisetype)];
    
    row = [ row ' fov: ' num2str(params.forces_on_v)];
    
    
end