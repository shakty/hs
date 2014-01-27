function wrapperanalysis(paramArgs)
    
    %% Add other directories to path
    path(path,'util/'); % Help functions
    path(path,'lib/'); % Help functions

    for i=1:length(paramArgs)
        params = (paramArgs{i});
        
        %truthradius_onefile(params);
        agents_onefile(params);
        clusters_onefile(params);
    end

end

