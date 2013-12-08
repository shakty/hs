function [] = temporalysis_truthradius_fun(dumpDir, subDir, simName)

    %% Saves simulations into properly formatted CSV files

    %% Add other directories to path
    path(path,'util/'); % Help functions
    path(path,'lib/'); % Help functions

    % Change default axes fonts.
    set(0,'DefaultAxesFontName', 'Times New Roman')
    set(0,'DefaultAxesFontSize', 14)

    % Change default text fonts.
    set(0,'DefaultTextFontname', 'Times New Roman')
    set(0,'DefaultTextFontSize', 14)


    % Dump to CSV only every X steps. 
    DUMP_RATE = 1; % 100 should produces ~22 iterations

    % Different radius to try out (around truth).
    RADIUSs = [0.01, 0.05, 0.1, 0.25, 0.4];
    
    % An agent must stay for at least STAY_FOR iterations.
    STAY_FOR = 20;
    
    % A share of at least CONSENSUS_THRESHOLD agents must stay within the
    % smallest radius to be say that there is a consensus on truth.
    CONSENSUS_THRESHOLD = 2/3;
    
    % Consenus must on truth must hold for at least CONSENSUS_ON_TRUTH_FOR
    % iterations to say that there is really a consensus on truth.
    CONSENSUS_ON_TRUTH_FOR = 20;
    
    CSV_DUMP = 1;
    PLOTS = 0;

    DUMPDIR = [dumpDir subDir];

    tic
    temporal_analysis_truthradius(DUMPDIR, simName, RADIUSs, STAY_FOR, ...
        CONSENSUS_ON_TRUTH_FOR, CONSENSUS_THRESHOLD, CSV_DUMP, ...
        DUMP_RATE, PLOTS);
    toc
end



