%fh = @truth_linear;
%fplot(fh,[0 1])

f = inline('x.*sin(x.*y)','x','y')

f = inline('x.*sin(x.*y)','x','y')

f = inline('y-(x- 0.5)','x','y')

[X,Y] = meshgrid(0:.01:1,0:20);
Z = f(X,Y)
%mesh(X,Y,Z)

surf(X,Y,Z)

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
score = v - A * (1 - exp( -d / d0))