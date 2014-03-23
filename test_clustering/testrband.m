clc
cla

BAND = 0.1; %0.09819;
vi = [0.1:0.05:0.5];
NELEM = length(vi);

vy = zeros(NELEM,1);

% for (i=1:NELEM)
%     vy(i) = rband(vi(i), BAND);
%     if (mod(i, 2) == 0)
%         color =  'r';
%     else
%         color = 'g';
%     end
%     cla
%     hold on
%     circle(0.5,0.5,vi(i),color);
%     circle(0.5,0.5,vy(i),color);
%     xlim([0,1]);
%     ylim([0,1]);
%     hold off    
%     pause(1)
% end

vi = 0;
vy = 0;

i = 1;

while (vi(i) <= 0.5)
    vy(i) = rband(vi(i), BAND, 0);
    if (mod(i, 2) == 0)
        color =  'r';
    else
        color = 'g';
    end
    %cla
    hold on
    circle(0.5,0.5,vi(i),color);
    circle(0.5,0.5,vy(i),color);
    xlim([0,1]);
    ylim([0,1]);
    hold off    
    %pause(1)
    i = i + 1;
    vi(i) = vy(i-1);
end

vi
% pause(10)
% 
% vi = 0.5;
% vy = 0;
% 
% i = 1;
% cla
% 
% while (vi(i) >= 0.05)
%     vy(i) = rband(vi(i), BAND, 1);
%     if (mod(i, 2) == 0)
%         color =  'r';
%     else
%         color = 'g';
%     end    
%     hold on
%     circle(0.5,0.5,vi(i),color);
%     circle(0.5,0.5,vy(i),color);
%     xlim([0,1]);
%     ylim([0,1]);
%     hold off    
%     %pause(1)
%     i = i + 1;
%     vi(i) = vy(i-1);
% end



% n = 1000;
% Xc = 0;
% Yc = 0;
% 
% theta = rand(1,n)*(2*pi);
% r = (rand(1,n)*(10-8)+8);
% x = Xc + r.*cos(theta);
% y = Yc + r.*sin(theta);
% plot(x,y,'.');
