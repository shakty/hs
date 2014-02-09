function [row] = create_params_string_small(params) 

    row = [ 'R: ' num2str(params.R) ' alpha: ' num2str(params.alpha) ];
        
    row = [ row ' tau: ' num2str(params.tau)];
    
    row = [ row ' sigma: ' num2str(params.sigma)];

    row = [ row ' epsilon: ' num2str(params.epsilon)];
    
    row = [ row ' v0: ' num2str(params.v_scaling)];

%     row = [ row ' agents: ' num2str(params.n_agents)];    
%     
%     row = [ row ' noise: ' num2str(params.noisetype)];
%     
%     row = [ row ' fov: ' num2str(params.forces_on_v)];
    
    
end