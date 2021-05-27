clc;
clear all;
close all;

f = 2.4e9;
c = 3e8;
lambda = c / f;
k = (2 * pi) / lambda;
psi = 0;
d = lambda / 2;
I = ones(5);
In = I(1, :);
theta = (0:0.1:180);

w = [1; 0.9018; 0.72759; 0.51502; 0.4159];
figure(1)
theta=[0:0.1:180];
L=length(theta);
for l=1:L
b=(2*(cos(pi/2*[1:2:9]*cos((pi/180)*theta(l)))));
AF(l)=abs(b*w);
end

N = 5;
I = ones(N);
In = I(1, :);
for(i = 1 : length(theta))
	%AF(:, i) = abs(2 * In * cos(k * d * cos(theta(i) * (pi / 180))));
	AF1(:,i)=abs(In*(2*(cos(pi/2*[1:2:9]*cos(theta(i)*(pi/180))))).');
end

plot(theta, 20 * log(AF / max(AF)));
hold on
plot(theta, 20 * log(AF1 / max(AF1)));
