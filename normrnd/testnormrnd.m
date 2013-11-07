clear
clc
N = 100000;
sigma = 0.5;

test = zeros(N,1);
v = 0;
for i=1:N
    test(i) = normrnd(0, sigma);
    v = v + test(i);
    
end
v
mean = mean(test)



v2 = [1;1];
norm(v2)
for i=1:N
    [THETA,RHO] = cart2pol(v2(1,1), v2(2,1));
    THETA2 = normrnd(THETA,sigma);
    RHO2 = normrnd(RHO,sigma);
    [X,Y] = pol2cart(THETA2,RHO2);
    v2(1,1) = X;
    v2(2,1) = Y;
    v2;
end

v2
norm(v2)


test = zeros(2,N);
v = [0;0];
for i=1:N
    test(:,i) = normrnd(v, sigma);
    v = test(:,i);
end
norm(v)