clear all
close all
clc
rng('shuffle')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1-Effect of Eigenvalue Spread
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 800;
K = 1000;
h1 = [0.2194, 1.0, 0.2194];
sigma_noise = 0.01;
M = 9;
step_size = 0.075;
e = run_ensemble(N, K, h1, sigma_noise, M, step_size, false);
figure('DefaultAxesFontSize',24);
semilogy(e, 'LineWidth', 2);
xlabel("Time Sample (n)")
ylabel("Normalized MSE")
hold on
h2 = [0.2798, 1.0, 0.2798];
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
semilogy(e, 'LineWidth', 2);
h3 = [0.3365 1.0 0.3365];
e = run_ensemble(N, K, h3, sigma_noise, M, step_size, false);
semilogy(e, 'LineWidth', 2);
legend('Channel 1', 'Channel 2', 'Channel 3','Channel 4');
title('Effect of Eigenvalue Spread')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, 'LMS Effect of Eigenvalue Spread.png');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2- Effect of Filter Order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2-i M = {9, 11, 21}, step = 0.075
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 9;
step_size = 0.075;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
figure('DefaultAxesFontSize',24);
semilogy(e, 'LineWidth', 2);
xlabel("Time Sample (n)")
ylabel("Normalized MSE")
hold on
M = 11;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
semilogy(e, 'LineWidth', 2);
M = 21;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
semilogy(e, 'LineWidth', 2);
legend('M = 9', 'M = 11', 'M = 21');
title('Effect of Filter Order at Step Size 0.075 (Channel 2)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, 'LMS Effect of Filter Order_i.png');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2-ii M = {9, 11, 21}, step = 0.0375
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 750;
M = 9;
step_size = 0.0375;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
figure('DefaultAxesFontSize',24);
semilogy(e, 'LineWidth', 2);
xlabel("Time Sample (n)")
ylabel("Normalized MSE")
hold on
M = 11;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
semilogy(e, 'LineWidth', 2);
M = 21;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
semilogy(e, 'LineWidth', 2);
legend('M = 9', 'M = 11', 'M = 21');
title('Effect of Filter Order at Step Size 0.0375 (Channel 2)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, 'LMS Effect of Filter Order_ii.png');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2-iii M = {9, 11, 21}, step = 0.0125
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 2500;
M = 9;
step_size = 0.0125;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
figure('DefaultAxesFontSize',24);
semilogy(e, 'LineWidth', 2);
xlabel("Time Sample (n)")
ylabel("Normalized MSE")
hold on
M = 11;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
semilogy(e, 'LineWidth', 2);
M = 21;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
semilogy(e, 'LineWidth', 2);
legend('M = 9', 'M = 11', 'M = 21');
title('Effect of Filter Order at Step Size 0.0125 (Channel 2)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, 'LMS Effect of Filter Order_iii.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3- Effect of Step-size Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1500; 
M = 9;
step_size = 0.0125;
e = run_ensemble(N, K, h1, sigma_noise, M, step_size, false);
figure('DefaultAxesFontSize',24);
semilogy(e, 'LineWidth', 2);
xlabel("Time Sample (n)")
ylabel("Normalized MSE")
hold on
step_size = 0.025;
e = run_ensemble(N, K, h1, sigma_noise, M, step_size, false);
semilogy(e, 'LineWidth', 2);
step_size = 0.075;
e = run_ensemble(N, K, h1, sigma_noise, M, step_size, false);
semilogy(e, 'LineWidth', 2);
legend('Step Size = 0.0125', 'Step Size = 0.025', 'Step Size = 0.075');
title('Effect of Step Size at Filter Order 9 (Channel 1)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, 'LMS Effect of Step Size.png');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4- Comparison of Standard LMS and Normalized LMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1500;
M = 9;
step_size = 0.025;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
figure('DefaultAxesFontSize',24);
semilogy(e, 'LineWidth', 2);
xlabel("Time Sample (n)")
ylabel("Normalized MSE")
hold on
step_size = 0.075;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, false);
semilogy(e, 'LineWidth', 2);
step_size = 1;
e = run_ensemble(N, K, h2, sigma_noise, M, step_size, true);
semilogy(e, 'LineWidth', 2);
legend('Step Size = 0.025', 'Step Size = 0.075', 'Normalized LMS');
title('Comparison of Standard LMS and Normalized LMS (Channel 2)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, 'LMS Comparison of Standard LMS and Normalized LMS.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function e = run_ensemble(N, K, h, sigma_noise, M, step, normalized_lms)
    e = zeros(K, N);
    for i = 1:K
        e(i,:) = run_experiment(N, h, sigma_noise, M, step, normalized_lms);
    end
    e = mean(e);
    %e = e / e(1);
end
 
function e = run_experiment(N, h, sigma_noise, M, step, normalized_lms)
    d=2*((randn(1,N)>0)-0.5)    
    x = conv(d, h);
    x = x(1:N);
    noise = sigma_noise*randn(1,N);
    u = x + noise;
    w = zeros(M, N);
    d_hat = zeros(1, N);
    e = zeros(1, N);
    delta1 = floor(0.5*length(h));
    delta2 = floor(0.5*M);
    delay = delta1 + delta2;
    for i = 1:N
        for j = 0:M-1
            if(i-j>0)
                d_hat(i) = d_hat(i) + w(j+1,i)*u(i-j);
            end
        end
        if (i>delay)
            e(i) = d(i-delay) - d_hat(i);
        end
        if(i<N)
            u_vector = zeros(M,1);
            u_vector(1:min(i, M)) = transpose(flip(u(max(1, i-M+1): i)));
            if(normalized_lms)
                step = 1/(norm(u_vector)).^2;
            end
            w(:,i+1) = w(:,i) + step*u_vector*e(i);
        end
    end
    e = e.^2;
end

