function [] = temporalysis_fun(dumpDir, subDir, simName)

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


    % Dump to CSV (micro and macro) only every X steps
    % the sums.mat is not affected by this value and always contains
    % all iterations
    DUMP_RATE = 100; % should produces ~22 iterations

    % Only clusters of size above the cutoff are included in the analysis
    CLU_CUTOFF = 2;
    % When computing the coverage we build a grid on top of the space of cell
    % size = PRECISION
    PRECISION = 100;

    CSV_DUMP = 1;
    PLOTS = 0;

    % DUMPDIR = 'dump/';
    % simName = 'test_t-2013-6-4-12-14/';

    DUMPDIR = [dumpDir subDir];

    tic
    temporal_analysis(DUMPDIR, simName, PRECISION, CLU_CUTOFF, CSV_DUMP, DUMP_RATE, PLOTS);
    toc
end



