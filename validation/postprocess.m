clear
clc
close all

data = csvread('../series/enstrophy.csv');
compdata = csvread('reference.csv');
figure
hold on
plot(data(:, 1), data(:, 2)/data(1, 2))
plot(compdata(:, 1), compdata(:, 2)/compdata(1, 2))
legend('me', 'them');
axis([0 20 0 14]);