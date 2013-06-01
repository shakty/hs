%fh = @truth_linear;
%fplot(fh,[0 1])

f = inline('x.*sin(x.*y)','x','y')

% tau = 0.2
% x position => x-0.5 = dist
% y speed
% if f < 0 the agent is attracted
f = inline('y-((x- 0.5)/0.2)','x','y')
[X,Y] = meshgrid(0:.01:1,0:20);
Z = f(X,Y)
surface(X,Y,Z)
mesh(X,Y,Z)

contour(X,Y,Z)
colorbar

f = inline('y-(x- 0.5)/0.2','x','y')
[X,Y] = meshgrid(0:.01:1,0:20);
Z = f(X,Y)
%surface(X,Y,Z)
%mesh(X,Y,Z)
contour(X,Y,Z)
colorbar



colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);
truth = [0.5;0.5]
tau = 0.1
[X,Y] = meshgrid(0:.01:1,0:.01:1);
a = [0:.01:1;0:.01:1];
n_agents = length(a)
ths = @(x,y) (repmat(truth,1,n_agents)-[x;y]).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-[x;y],2),2,1)).*normpdf(abs(repmat(truth,1,n_agents)-x),norm(truth)./2,0.1);
for i=1:length(X)
    Z(i) = ths(X(i),Y(i))
end
contour(X,Y,Z)
colorbar


ths = @(x,y) (repmat(truth,1,n_agents)-[x;y])

y = (repmat(truth,1,n_agents)-x)./(norm((repmat(truth,1,n_agents)-x))./tau);

clear all; close all;
MIN = realmin('single')
syms x 
y= (1/log(x) + 14) / 14% - log(MIN) - 83.5;
h=ezplot(y,[0,1]);
set(h, 'Color', 'r'); 
hold on;
y=abs(x);
h=ezplot(y,[0,1]);
set(h, 'Color', 'r'); 
%hold on;
%ezplot(diff(y,x),[0,1]);


d = 0.1;
v = 0;
d0 = 1;
A = 1;
score = v - A * (1 - exp( -d / d0));

colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);
a = [0:0.05:1;0:0.05:1]
truth = [0.5;0.5];
n_agents = length(a);
tau = 0.1
% NORMAL
ths = @(x) (repmat(truth,1,n_agents)-x).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-x,2),2,1)).*normpdf(abs(repmat(truth,1,n_agents)-x),norm(truth)./4,0.1);
ths(a)

% TRUTH Accelerating closer to Truth (LOG-NORMAL)
ths = @(x) (repmat(truth,1,n_agents)-x).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-x,2),2,1)).*lognpdf(abs(repmat(truth,1,n_agents)-x),log(norm(truth)),1);
ths(a)




x = (10:1000:125010)';
y = lognpdf(x,log(20000),1.0);
plot(x,y)
set(gca,'xtick',[0 30000 60000 90000 120000])
set(gca,'xticklabel',{'0','$30,000','$60,000',...
                             '$90,000','$120,000'})
                         
                         
mean(y)