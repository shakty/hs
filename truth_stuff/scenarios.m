clc;
clear;

path(path,'../util/'); % Help functions

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 14)

DIAG = sqrt(2); 
colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);

PRECISION = 0.01;
x = [0:PRECISION:1;0:PRECISION:1];


tau = 0.1;

truth = [0.5;0.5];
truth = [0.1;0.1];


% TRUTH Constant
ths = @(x) (repmat(truth,1,length(x))-x)./tau.*(repmat(tau./colnorm(repmat(truth,1,length(x))-x,2),2,1));

% TRUTH Linear
ths = @(x) (repmat(truth,1,length(x))-x)./tau;

% TRUTH Decaying exponentially (EXP)
SIGMA = 1;
ths = @(x) (repmat(truth,1,length(x))-x)./tau.*repmat(exppdf(colnorm(repmat(truth,1,length(x)) - abs((repmat(truth,1,length(x))-x)),2),SIGMA),2,1);

% Millean Arena (NORMAL)
POS = 3; SIGMA = 0.05;
ths = @(x) (repmat(truth,1,length(x))-x)./tau.*repmat(normpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);

% Hard to Find (NORMAL)
POS = 100; SIGMA = 0.02;
ths = @(x) (repmat(truth,1,length(x))-x)./tau.*repmat(normpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);

% Wide Funnel to Truth (LOG-NORMAL)
POS = 1; SIGMA = 3;
ths = @(x) (repmat(truth,1,length(x))-x)./tau.*repmat(lognpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);
 
% Gentle Landing to truth
POS = 0; SIGMA = 0.2;
ths = @(x) (repmat(truth,1,length(x))-x)./tau.*repmat(normpdf(colnorm(repmat(truth,1,length(x))-x,2),0,0.2),2,1);
        

str = 'Wide Funnel to Truth';
tp = '_tm';
h = figure;
[X,Y] = meshgrid(x(1,:), x(1,:));
Z = zeros(size(X));
for i=1:length(X)
    agents = [X(i,:) ; Y(i,:)];
    forces = ths(agents);
    Z(i,:) = colnorm(forces,2);
end
subplot(2,2,1:2);
mesh(X,Y,Z);
%colorbar;
title('3D');
subplot(2,2,3);
contour(X,Y,Z);
hold on
plot(truth(1), truth(2),'rx');
hold off
title('2D');
subplot(2,2,4);
hold on
line = ths(x);
plot(x, line(1,:), 'b');
hold on
plot(truth(1), 0,'ro');
hold off
title('1D');
%Adding suptitle
suptitle(str);
str = strrep(str, ' ', '_');
saveas(h,[ 'scenarios/' str tp '.fig']);
saveas(h,[ 'scenarios/' str tp '.png']);


hold on
cdfplot(abs(line(1,:)));
hold off

title(['pos: ' num2str(POS) ', sigma: ' num2str(SIGMA)]);


prompt = 'Give a file name to save the image.';
str = input(prompt,'s');
if ~isempty(str)
    saveas(h,[ 'scenarios/' str ]);
end

return;


%%% To find parameters;

for SIGMA=0.1:0.1:20

    for POS=0.1:0.1:8
       
        ths = @(x) (repmat(truth,1,length(x))-x)./tau.*repmat(normpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);
        
        [X,Y] = meshgrid(x(1,:), x(1,:));
        Z = zeros(size(X));
        for i=1:length(X)
            agents = [X(i,:) ; Y(i,:)];
            forces = ths(agents);
            Z(i,:) = colnorm(forces,2);
        end
        h = mesh(X,Y,Z);
        
        hold on
        plot(truth(1), truth(2),'rx');
        hold off
        title(['pos: ' num2str(POS) ', sigma: ' num2str(SIGMA)]);
       
       
        prompt = 'Give a file name to save the image.';
        str = input(prompt,'s');
        if ~isempty(str)
            saveas(h,[ 'scenarios/' str ]);
        end
        
    end
end
return;

contour(X,Y,Z);
colorbar;
surface(X,Y,Z);



return;

hold on
%plot(x, const(1,:), 'r');
%plot(x, linear(1,:), 'k-')
plot(x, normMid(1,:))
%plot(x, expo(1,:))
%plot(x, 15*normMoved(1,:))
%plot(x, lognorm(1,:));
%plot(x, 20*expo(1,:));
%plot(x, lognorm(1,:));
line(truth,[-1, 1]);
line([0, 1],[0,0]);
hold off

hold on
%cdfplot(abs(const(1,:)));
%cdfplot(abs(linear(1,:)));
cdfplot(abs(normMid(1,:)));
%cdfplot(abs(normMoved(1,:)));
%cdfplot(abs(lognorm(1,:)));
hold off
