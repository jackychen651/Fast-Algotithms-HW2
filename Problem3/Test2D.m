clear
clc
close all
% We assume that all x will be in (-3,-1)*(-1,1) and y will be in (1,3)*(-1,1)
% Also q_x is the charges of x, q_y is the charges of y.
% M is truncation terms
% N is number of charges in one box
rng(42);

N=1000;
M=30;
x_real=rand(N,1)*2+3;
x_im=rand(N,1)*2-1;
[theta_x, r_x] = cart2pol(x_real, x_im);
y_real=rand(N,1)*2-1;
y_im=rand(N,1)*2-1;
[theta_y, r_y] = cart2pol(y_real, y_im);
q_y=rand(N,1)*2-1;
u = FMM2D(r_x, theta_x, r_y, theta_y, q_y, M);
u_truth = GroundTruth2D(r_x, theta_x, r_y, theta_y, q_y);
error = norm(u-u_truth,2)/sum(abs(q_y),1);

%% test
t1=[];
t2=[];
errors=[];
Ns=[];
N=32;
P=floor(log(eps)/log(sqrt(2)/3));
i=1;
while N < 10000
    x_real=rand(N,1)*2+3;
    x_im=rand(N,1)*2-1;
    [theta_x, r_x] = cart2pol(x_real, x_im);
    y_real=rand(N,1)*2-1;
    y_im=rand(N,1)*2-1;
    [theta_y, r_y] = cart2pol(y_real, y_im);
    q_y=rand(N,1)*2-1;

    tic
    u = FMM2D(r_x, theta_x, r_y, theta_y, q_y, M);
    run_time1=toc;
    t1=[t1;run_time1];

    tic
    u_truth = GroundTruth2D(r_x, theta_x, r_y, theta_y, q_y);
    run_time2=toc;
    t2=[t2;run_time2];

    error = norm(u-u_truth,2)/sum(abs(q_y),1);
    errors = [errors,error];
    
    Ns=[Ns,N];
    fprintf("i=%d,N=%d,t1=%e,t2=%e,error=%e\n", i,N,t1(i),t2(i),errors(i));
    %fprintf("i=%d,P=%d,t1=%e\n", i,P,t1(i));
    N=N*2;
    i=i+1;
    save('Ns2D.mat', 'Ns');
    save('t12D.mat', 't1');
    save('t22D.mat', 't2');
    save('errors2D.mat', 'errors');
end

%%
Ns = load('Ns2D.mat');
t1 = load('t12D.mat');
t2 = load('t22D.mat');
errors = load('errors2D.mat');
Ns=Ns.Ns;
t1=t1.t1;
t2=t2.t2;
errors=errors.errors;

figure;
plot(Ns, t1, '-o', 'DisplayName', 'Fast Multipole Method');
hold on;
plot(Ns, t2, '-s', 'DisplayName', 'Direct Method');
hold off;
xlabel('Number of Charges in One Box');
ylabel('Time');

legend('show', 'Location', 'northwest')

grid on;
set(gca, 'LooseInset', get(gca, 'TightInset'));
width = 6;
height = 6;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperSize', [width height]);
print(gcf, 'Benchmark2D', '-dpdf', '-fillpage');

%%
figure;
plot(Ns, errors, '-o', 'DisplayName', 'Relative Error');
hold on;

xlabel('Number of Charges in One Box');
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
print(gcf, 'Error2D', '-dpdf', '-fillpage');

