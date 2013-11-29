function wrappersim( paramArgs )
    
    %% Add other directories to path
    path(path,'util/'); % Help functions
    path(path,'lib/'); % Help functions

    for i=1:length(paramArgs)
        simulation(paramArgs{i})
    end

end

