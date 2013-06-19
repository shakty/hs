%% Saves simulations into properly formatted CSV files

close all;
clear;
clc;

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 14)

% The minimum size of a cluster to be included in the analysis
CLU_CUTOFF = 10;
% When computing the coverage we build a grid on top of the space of cell
% size = PRECISION
PRECISION = 100;

CSV_DUMP = 1;
PLOTS = 0;

DUMPDIR = 'dump/';

<<<<<<< HEAD
simName = 'attrExpo_nv_rndseq_tm_Rleft/'; 

=======
simName = 'test_t-2013-6-4-12-14/';
 
dumpDir = [DUMPDIR simName '/'];
tic  
files = dir(dumpDir);
fileIndex = find(~[files.isdir]);

validFileIdx = 0;
for f = 1:length(fileIndex)

        append = files(fileIndex(f)).name;
        fileName = [dumpDir, append];

        % We load only .mat
        [PATH,NAME,EXT] = fileparts(fileName);
        if (~strcmpi(EXT,'.mat') || strcmp(NAME, 'temporalysis') == 1) 
            continue;
        end
        
        simnameidx = strfind(NAME, '-');
        simnameidx = NAME(1:simnameidx-1);
        
        
        load(fileName);
        validFileIdx = validFileIdx + 1
end
toc
return;
>>>>>>> 77e9dabdfd177e769922723a616214ac001cb06c
%profile on
tic
temporal_analysis(DUMPDIR, simName, PRECISION, CLU_CUTOFF, CSV_DUMP, PLOTS);
toc
%profile viewer
%%




