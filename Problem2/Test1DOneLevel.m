clear;
clc;
close all;

profile on
rng(42);
% use customary random seed

% Notations and Assumptions
% u is interaction, x is point coordinates, q is charges
% Assume we have P points, and we divide the intervals to M parts
% We also assume that all points are within [-1, 1], unifrom distribution
% And the charges are also with [-1, 1], unifrom distribution

% P = 2;
% x = -1+2*rand(P, 1);
% q = -1+2*rand(P, 1);
% k = 10;
% u_true = GroundTruth1D(x,q,k);

%% test
t1=[];
t2=[];
errors=[];
Ps=[];
K=10;
P=32;
i=1;
while P < 100000  %100000
    x = -1+2*rand(P, 1);
    q = -1+2*rand(P, 1);
    tic
    u = FMM1DOneLevel(x,q,K);
    run_time1=toc;
    t1=[t1;run_time1];

    tic
    u_true = GroundTruth1D(x,q,K);
    run_time2=toc;
    t2=[t2;run_time2];

    error = norm(u - u_true, 2) / sum(abs(q),1);
    errors = [errors,error];
    
    Ps=[Ps;P];
    fprintf("i=%d,P=%d,t1=%e,t2=%e,error=%e\n", i,P,t1(i),t2(i),errors(i));
    P=P*2;
    i=i+1;
end
save('PsOneLevel.mat', 'Ps');
save('t1OneLevel.mat', 't1');
save('t2OneLevel.mat', 't2');
save('errorsOneLevel.mat', 'errors');

%% Plot

Ps = load('PsOneLevel.mat');
t1 = load('t1OneLevel.mat');
t2 = load('t2OneLevel.mat');
errors = load('errorsOneLevel.mat');
Ps=Ps.Ps;
t1=t1.t1;
t2=t2.t2;
errors=errors.errors;

figure;
plot(Ps, t1, '-o', 'DisplayName', 'Fast Multipole Method');
hold on;
plot(Ps, t2, '-s', 'DisplayName', 'Direct Method');
hold off;

xlabel('Number of Charges');
ylabel('Time');

legend('show', 'Location', 'northwest')

grid on;
set(gca, 'LooseInset', get(gca, 'TightInset'));
width = 6;
height = 6;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperSize', [width height]);
print(gcf, 'Benchmark1DOneLevel', '-dpdf', '-fillpage');


%%
figure;
plot(Ps, errors, '-o', 'DisplayName', 'Relative Error');
hold on;

xlabel('Number of Charges');
ylabel('Relative Error');

set(gca, 'YScale', 'log');

legend('show', 'Location', 'northwest')

grid on;
set(gca, 'LooseInset', get(gca, 'TightInset'));
width = 6;
height = 6;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperSize', [width height]);
print(gcf, 'Error1DOneLevel', '-dpdf', '-fillpage');