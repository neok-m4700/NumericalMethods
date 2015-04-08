npoints = 200;
%linearly increasing x
xs = linspace(1,100,npoints)';
m = rand*20 - 20;
c = rand*100 - 50;
% linearly increasing y
ys_clean = m*xs + c;

noise_factor = xs+100; %scaled noise
noise = rand(npoints,1).*noise_factor - noise_factor/2;
%noisy y
ys = ys_clean + noise;

%solve for m and c
points = [xs ys];
[m c] = linfit(points);
ys_solved = m*xs + c;

%plot
scatter(xs,ys,'r') 
hold on
plot(xs,ys_solved);
hold off
