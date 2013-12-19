clear
clc
N = 100000;
sigma = 0.5;

colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);


% test = zeros(N,1);
% test2 = zeros(N,1);
% v = 0;
% v2 = 0;
% for i=1:N
%     test(i) = normrnd(0, sigma);
%     v = v + test(i);
%     if (i>1)
%         test2(i) = normrnd(test2(i-1), sigma);
%     end
%     v2 = v2 + test2(i);
% end
% v
% mean1 = mean(test)
% 
% v2
% mean2 = mean(test2)

% test = zeros(N,1,2);
% test2 = zeros(N,1,2);
% v = 0;
% v2 = 0;
% 
% for i=1:N
%     test(i,:) = normrnd([0;0], sigma, 2, 1);
%     v = v + test(i);
%     if (i>1)
%         test2(i,:) = normrnd(test2(i-1), sigma, 2, 1);
%     end
%     v2 = v2 + norm(test2(i,:));
% end
% v
% mean1 = mean(test)
% 
% v2
% mean2 = mean(test2)




v2 = [1;1];
norm(v2)
for i=1:N
    [THETA,RHO] = cart2pol(v2(1,1), v2(2,1));
    THETA2 = normrnd(THETA,sigma);
    RHO2 = normrnd(RHO,sigma);
    [X,Y] = pol2cart(THETA2, RHO);
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