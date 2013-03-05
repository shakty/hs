%% Produces the HEATMAP

close all;
clear;
clc;

CONV_PLOT=1;    %select one parameter
HEAT_MAP=1;     %select two parameters
MOVIE=1;

DUMPDIR = 'dump/';

%Data Format: conv,t_end,dt,agents,iss,isd,A,B,k,d0,d1,alpha,tau,R,sigma,v_scaling,nof_cluster,clusterTightness;
nPARAMETERS = 18; % 17+conv

%% Load Data

%dumpDir = [DUMPDIR '2011-5-24-11-24/'];
dumpDir = [DUMPDIR '2011-5-24-11-57/'];
%dumpDir = [DUMPDIR '2011-5-24-15-4/'];
%dumpDir = [DUMPDIR '2011-5-24-19-31/'];

simdir = 'few_big_groups-DIM-vs-ALPHA';

dumpDir = [DUMPDIR simdir '/'];




files = dir(dumpDir);
fileIndex = find(~[files.isdir]);

data = zeros(size(fileIndex,1),nPARAMETERS);

if (isempty(fileIndex))
    error('Invalid Directory Selected');
end

% Load all parameters matrices in one
for i = 1:length(fileIndex)
    append = files(fileIndex(i)).name;
    fileName = [dumpDir, append];
    load(fileName);
    %data(i,:) = [dump{4}, dump{1}];
    data(i,:) = [cell2mat(struct2cell(dump.parameters)); dump.conv];
end
% All convergences are saved in one array;
convs = data(:,1);



%% Select Data
%Which parameters shall be fixed and at what value? 
%negative value means we do not care,
%i.e. the parameter is irrelevant (-2) or to be under examination (-1)

conv_mes = -2;  
t_end   = -2;
dt      = -2;
agents  = 100;
iss     = -2;
isd     = -2;
A       = -1;
B       = -1;
k       = 1;
d0      = 1;
d1      = 1;
alpha   = 1;
tau     = 1;%-1;
R       = 2;
sigma   = 0;
v_scaling = 0;
nof_cluster = 3;
clusterTightness = -2;

parameters = [conv_mes,t_end,dt,agents,iss,isd,A,B,k,d0,d1,alpha,tau,R,sigma,v_scaling,nof_cluster,clusterTightness];
param_names = {'Convergence Measure','t_end','dt','nof. agents','space size',...
    'space dim.','A','B','k','d0','d1','alpha','tau','R','sigma',...
    'init. velocities scaling','nof. init. clusters','init. cluster radii scaling'};

selected_parameters = [];
fixed_parameters = [];
eval_data = [];


for i=1:nPARAMETERS
    if(parameters(i)>=0) %the parameter is restricted: select relevant rows
        fixed_parameters = [fixed_parameters,i];
    else if(parameters(i)==-1) %the parameter shall be examined 
            selected_parameters = [selected_parameters,i];
            eval_data = [eval_data, data(:,i)];
        end
    end
end


%% Convergence over Parameter Plot
if(CONV_PLOT)
    if(size(selected_parameters,2)==1)
        figure;
        plot(eval_data(:,selected_parameters),eval_data(:,1),'k-x','LineWidth',2);
        set(gca,'FontSize',15);
        fixed_string = ['Fixed Param.:'];
        for i=fixed_parameters(2:end-2)
            fixed_string = [fixed_string ' ' param_names{i} '=' int2str(parameters(i)) ];
        end
        title_string = sprintf(['Convergence in dependance of ' param_names{selected_parameters} '\n' fixed_string]);
        title(title_string);
        xlabel(['Parameter ' param_names{selected_parameters}]);
        ylabel('Quadratic Mean of Dist. to Truth');
        ylim([0,0.2]);
        grid on;
        box on;
    end
end

%% Phase Diagram (Heat Map)
if(HEAT_MAP)
     if(size(selected_parameters,2)==2)
        figure;
        set(gca,'FontSize',15);    
        
        combinations = unique(eval_data,'rows');
        
        x=eval_data(:,1);
        y=eval_data(:,2);
        
        xu=unique(x);
        yu=unique(y);
        
        [a,b]=meshgrid(xu,yu);
        
        M=zeros(size(xu,2),size(yu,2));
        
        
        for i=1:size(eval_data,1)
            find(eval_data == combinations);
            
            v1=eval_data(i,1);
            v2=eval_data(i,2);
            M(find(xu==v1),find(yu==v2))=convs(i);
        end
        
        M;
        
        hold on;
        view(0,90);
        surf(a',b',M,'EdgeColor','none');
                
        xlabel(['Parameter ' param_names{selected_parameters(1)}]);
        ylabel(['Parameter ' param_names{selected_parameters(2)}]);
        fixed_string = 'Fixed Param.:';
        
        for i=fixed_parameters(2:end-2)
            fixed_string = [fixed_string '  ' param_names{i} '=' int2str(parameters(i)) ];
        end
        
        title_string = sprintf(['Phase Diagram: Convergence in Dependance of '...
            param_names{selected_parameters(1)} ' and ' param_names{selected_parameters(2)}...
            '\n' fixed_string]);
        title(title_string);       
        
        colorbar;
        box on;
     end
end

%% Movie Rendering
% if(MOVIE)
%     mov = avifile('movie.avi');
%     for i=1:100
%         plot();
%         F = getframe(gca);
%         mov = addframe(mov,F);
%     end
%     mov = close(mov);
% end
