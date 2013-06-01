path(path,'../util/'); % Help functions

d = 0.1;
v = 0;
d0 = 1;
A = 1;
score = v - A * (1 - exp( -d / d0));

colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);
PRECISION = 0.001;
a = [0:PRECISION:1;0:PRECISION:1];
truth = [0.5;0.5];
n_agents = length(a);
tau = 0.1;

% TRUTH Constant
ths = @(x) (repmat(truth,1,n_agents)-x).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-x,2),2,1));
const = ths(a);

% TRUTH Linear
ths = @(x) (repmat(truth,1,n_agents)-x)./tau;
linear = ths(a);

% TRUTH Accelerating in the middle (NORMAL)
ths = @(x) (repmat(truth,1,n_agents)-x).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-x,2),2,1)).*normpdf(abs(repmat(truth,1,n_agents)-x),norm(truth)./2,0.1);
normMid = ths(a);

% TRUTH Accelerating in the middle (NORMAL) mean closer to TRUTH
ths = @(x) (repmat(truth,1,n_agents)-x).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-x,2),2,1)).*normpdf(abs(repmat(truth,1,n_agents)-x),norm(truth)./4,0.1);
normMoved = ths(a);

% TRUTH Accelerating closer to Truth (LOG-NORMAL)
ths = @(x) (repmat(truth,1,n_agents)-x).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-x,2),2,1)).*lognpdf(abs(repmat(truth,1,n_agents)-x),log(norm(truth)),1);
lognorm = ths(a);

% TRUTH Decaying exponentially (EXP)
%ths = @(x) exppdf(abs(repmat(truth,1,n_agents) - abs((repmat(truth,1,n_agents)-x)))./tau);
ths = @(x) (repmat(truth,1,n_agents)-x).*(repmat(tau./colnorm(repmat(truth,1,n_agents)-x,2),2,1)).*exppdf(abs(repmat(truth,1,n_agents) - abs((repmat(truth,1,n_agents)-x)))./tau);
expo = ths(a);

% TRUTH Decaying logaritmically (EXP)
ths = @(x) (repmat(truth,1,n_agents)-x).*exppdf(abs(repmat(truth,1,n_agents)-x),norm(truth));
loggo = ths(a);


% TRUTH Decaying logaritmically (EXP)
ths = @(x) (repmat(truth,1,n_agents)-x).*lognpdf(abs(repmat(truth,1,n_agents)-x),norm(truth));
loggo = ths(a);


% TRUTH Decaying logaritmically (EXP)
ths = @(x) (repmat(truth,1,n_agents)-x).*normpdf(abs(repmat(truth,1,n_agents)-x),norm(truth)./4,0.1);
loggo = ths(a);

% NORMAL
ths = @(x) (repmat(truth,1,n_agents)-x).*normpdf(abs(repmat(truth,1,n_agents)-x),norm(truth)./2,0.1);
normMid = ths(a);

% NORMAL MOVED
ths = @(x) (repmat(truth,1,n_agents)-x).*normpdf(abs(repmat(truth,1,n_agents)-x),norm(truth)./4,0.1);
normMoved = ths(a);


% TRUTH Accelerating closer to Truth (LOG-NORMAL)
ths = @(x) (repmat(truth,1,n_agents)-x).*lognpdf(abs(repmat(truth,1,n_agents)-x),log(norm(truth)),1);
lognorm = ths(a);

lognorm = ths(a);

X = a(1,:);
Y = a(:,1);
Z = ths(a);

surf(X,Y,Z);

x = reshape(X,length(X),length(X))';

y = reshape(Y,5,5)';
z = reshape(Z,5,5)';



f = inline('y-((x- 0.5)/0.2)','x','y')
[X,Y] = meshgrid(0:.01:1,0:20);
Z = f(X,Y)
surface(X,Y,Z)
mesh(X,Y,Z)

y = 1
x = 1
truth - [x;y]


ths = @(x,y) (repmat(truth,length(x))-[x,y])./tau;

ths = @(x,y) x-y;

aa = ths(X,Y);

[X,Y] = meshgrid(a);





for i=1:size(X,2)
    Z
end

surface(X,Y,Z)
mesh(X,Y,Z)


[D P] = allfitdist(10*lognorm(1,:),'PDF')



D(1)
aa = normpdf(a,0.5,0.2)
plot(a,aa)
plot(a, log(aa));

figure

hold on
%plot(a, const(1,:), 'r');
%plot(a, linear(1,:), 'k-')
%plot(a, 15*normMid(1,:))
%plot(a, 15*normMoved(1,:))
%plot(a, lognorm(1,:));
%plot(a, 20*expo(1,:));
plot(a, lognorm(1,:));
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

cdf = a;
for i = 1:length(a)
    cumSum = sum(normMid(1,1:i));
    if isnan(cumSum) 
        cdf(i) = cdf(i-1);
    else
        cdf(i) = cumSum;
    end
end

plot(a, cdf(1,:))

ecdf(abs(normMid(1,:)));

sum(abs(linear(1,:)))

%legend(h, 'logNorm', 'normMoved', 'normMid', 'const');



x = (0:0.02:10);
y = lognpdf(x,0,0.5);
plot(x,y); grid;
xlabel('x'); ylabel('p')

LIMIT=1000
ff = lognpdf(a, 0, 1);
plot(a(1:LIMIT), ff(1:LIMIT));

y = fpdf(a,3,1);
plot(a,y(1,:)/(0.1))