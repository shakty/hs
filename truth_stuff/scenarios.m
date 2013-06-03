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

truth = [0.5;0.5];
truth = [0.1;0.1];

tau = 0.1;

% TRUTH Constant
ths = @(x) (repmat(truth,1,length(x))-x).*(repmat(tau./colnorm(repmat(truth,1,length(x))-x,2),2,1));
const = ths(x);

% TRUTH Linear
ths = @(x) (repmat(truth,1,length(x))-x)./tau;
linear = ths(x);

% TRUTH Decaying exponentially (EXP)
ths = @(x) (repmat(truth,1,length(x))-x).*repmat(exppdf(colnorm(repmat(truth,1,length(x)) - abs((repmat(truth,1,length(x))-x)),2),SIGMA),2,1);
expo = ths(x);

% TRUTH Accelerating in the middle (NORMAL)
POS = 6; SIGMA = 0.4;
ths = @(x) (repmat(truth,1,length(x))-x).*repmat(normpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);
normMoved = ths(x);

% Accelerating close to the TRUTH
POS = 3; SIGMA = 2;
ths = @(x) (repmat(truth,1,length(x))-x).*repmat(normpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);
normMoved = ths(x);

% SUPER-THIN Funnel to the Truth (LOG-NORMAL)
POS = 1; SIGMA = 3;
ths = @(x) (repmat(truth,1,length(x))-x).*repmat(lognpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);
lognorm = ths(x);
 
% THIN VULCANO (does it work with truth in the corner??)
ths = @(x) (repmat(truth,1,length(x))-x).*repmat(normpdf(colnorm(repmat(truth,1,length(x))-x,2),0.01,0.1),2,1);
normMid = ths(x);

% LOGNORMAL
ths = @(x) (repmat(truth,1,length(x))-x).*repmat(lognpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);
        

for SIGMA=0.1:0.1:20

    for POS=6:20
       
        %POS = 6; SIGMA = 0.3;
        ths = @(x) (repmat(truth,1,length(x))-x).*repmat(normpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);

         
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
