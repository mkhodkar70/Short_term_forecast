clc
clear all
close all

%% Variables

F = 8;  %% forcing

r = 700;           %% reduced order of the main data (u)
q = 20;            %% embedded delay for Hankel representation
J = 2;             %% # of neighboring points
p = (J+1)*r;       %% reduced order of the control (u^2)

text = ['train/Lorenz96_F_', num2str(F), '.mat'];
load(text, 'X', 'tau');

X_train = X(:, 1:100000);

clear X

control_train = nonlinear_observables(X_train, J);

nts = size(X_train, 2);

%% Training 

[A_dmdc, B_dmdc, phi_dmdc, Lambda_dmdc, U_dmdc, U_1, U_2] = HDMDC_Execution(X_train, control_train, p, r, q, tau);

Lambda_dmdc = diag(Lambda_dmdc);

cnt_blow = 0;    %% number of eigenvalues out of unit circle

for(i = 1:r)
    if(abs(Lambda_dmdc(i)) > 1)
        cnt_blow = cnt_blow + 1;
    end
end

Lambda_dmdc = diag(Lambda_dmdc);

text = ['DMDC_F_', num2str(F), '_r_', num2str(r), '_p_', num2str(p), '_q_', num2str(q), '.mat'];
save(text, 'A_dmdc', 'B_dmdc', 'U_dmdc', 'U_1', 'U_2', 'phi_dmdc', 'Lambda_dmdc', 'phi_dmdc', 'nts', 'cnt_blow');

%% Variables

if(F == 4)
    Lambda_max = 0.54;
    D_KY       = 15;
elseif(F == 8)
    Lambda_max = 1.72;
    D_KY       = 28;
elseif(F == 16)
    Lambda_max = 2.55;
    D_KY       = 32;
end

m_test = 2000;                 %% length of testing set
m_i    = 22000;                %% starting point of testing set
m_f    = m_i + m_test - 1;     %% final point of testing set

load(['test/Lorenz96_F_', num2str(F), '.mat'], 'X', 'tau');

n_test = 10;
sd     = sqrt(var(X(n_test, :)));

t_i = 0;                       %% initial time of prediction      
t_f = t_i + (m_test-1)*tau;    %% final time of prediction

X_test = X(:, m_i:m_f);

m_test = size(X_test, 2);
n      = size(X_test, 1); 

x = 1:n;

t  = t_i:tau:t_f;        %% time
t  = Lambda_max*t;

%% Prediction via Hankel-DMD to increase the length of dataset

load(['DMDC_F_', num2str(F), '_r_', num2str(r), '_p_', num2str(p), '_q_', num2str(q), '.mat']);

X_pred = [];
x_0    = X_test(:, 1:q);
X      = Delay(x_0, n, q);
x_0_r  = regress(X, U_dmdc);

control = nonlinear_observables(x_0, J);
control = Delay(control, size(control, 1), q);

for(j = 1:m_test-q)
    
    if(mod(j, 10) == 0)
        j
    end
           
    x_new  = A_dmdc*x_0_r + B_dmdc*control;
    
    x_0_r = x_new;
    
    X_new        = U_dmdc*x_new;
    X_pred(:, j) = X_new(1:n); 
    
    X_temp = zeros(n, q);
    for(i = 1:q)
        X_temp(:, i) = X_new(n*(i-1)+1:i*n);         %% unfolding the delay-embedded vector, i.e. a nq column vector will be arranged in the form of a n*q matrix
    end
    control      = nonlinear_observables(X_temp, J);
    control      = Delay(control, size(control, 1), q); 
end

t = t(:, q+1:end);
t = t - t(1);

clc

%% Error Calculation

E = zeros(size(X_pred, 2), 1);

for(i = 1:length(E))
    E(i) = norm(X_pred(:, i) - X_test(:, i))/norm(X_test(:, i));
end

i_blow = floor(m_test/1.5)-q;
for(i = 1:length(E))
    if(E(i) > 0.3)
        i_blow = i;
        break
    end
end

t_blow = t(i_blow)
E_ave  = mean(E(1:i_blow))

%% Plots

figure();
plot(t(1:ceil(1.5*i_blow)), X_pred(n_test, 1:ceil(1.5*i_blow))/sd, '--', 'linewidth', 2, 'color', 'red');
hold on
plot(t(1:ceil(1.5*i_blow)), X_test(n_test, 1:ceil(1.5*i_blow))/sd, 'linewidth', 2, 'color', 'blue');
hold on
plot(t_blow*ones(100, 1), linspace(-4, 4, 100), ':', 'linewidth', 2, 'color', 'black');
xlim([0 12]);
ylim([-3 3])
hold on
legend('Prediction', 'Truth');
text = 'Time series for anomalous variable at $j$ = %d'; 
text = sprintf(text, x(n_test));
title(text, 'interpreter', 'latex');
xlabel('$\Lambda_{max} t$', 'interpreter', 'latex', 'Fontsize', 14);
ylabel('$X_{18}/\sigma$', 'interpreter', 'latex', 'Fontsize', 14);
set(gca, 'Fontsize', 13);


figure();
plot(t, 100*E, 'linewidth', 2, 'color', 'red');
hold on
plot(t_blow*ones(100, 1), linspace(0, 120, 100), ':', 'linewidth', 2, 'color', 'black');
hold on
xlim([0 12]);
ylim([0 120]);
title('Relative Error');
xlabel('$\Lambda_{max} t$', 'interpreter', 'latex', 'Fontsize', 14);
ylabel('$E(\%)$', 'interpreter', 'latex', 'Fontsize', 14);
set(gca, 'Fontsize', 13);



i_test = [75 158 252 440 670 790];
t_test = t(i_test);

figure();
for(j = 1:length(t_test))
    subplot(2, 3, j)
    plot(x, X_test(:, i_test(j)), 'linewidth', 2, 'color', 'blue');
    hold on
    plot(x, X_pred(:, i_test(j)), '--', 'linewidth', 2, 'color', 'red');
    hold on
    pbaspect([2, 1, 1]);
    text = '%.2f';
    text = sprintf(text, t_test(j));
    xlim([1 n]);
    title(['$\Lambda_{max} t =$ ', text]);
    xlabel('$j$', 'interpreter', 'latex', 'Fontsize', 14);
    ylabel('$X$', 'interpreter', 'latex', 'Fontsize', 14);
    set(gca, 'Fontsize', 14);
end



figure();

t_contour = 1200;

subplot(3, 1, 1);
show = X_test(:, 1:t_contour);
[CR, h] = contourf(t(1:t_contour), x, show, 'linestyle', 'none');
hold on
colormap(b2r(-12, 12));
cb1 = colorbar;
caxis([-12 12]);
pbaspect([4, 1, 1]);
xlim([0 16]);
ylim([1 n]);
% set(gca, 'YTick', [0:4:16])
xlabel('$\Lambda_{max} t$', 'interpreter', 'latex');
ylabel('$j$', 'interpreter', 'latex');
title('$\mathrm{Truth}$', 'interpreter', 'latex', 'Fontsize', 14);
set(gca, 'Fontsize', 12)


subplot(3, 1, 2);
show = X_pred(:, 1:t_contour);
[CR, h] = contourf(t(1:t_contour), x, show, 'linestyle', 'none');
hold on
plot(t_blow*ones(100, 1), linspace(1, n, 100), ':', 'linewidth', 2, 'color', 'black');
hold on
colormap(b2r(-12, 12));
cb2 = colorbar;
caxis([-12 12]);
pbaspect([4, 1, 1]);
xlim([0 16]);
ylim([1 n]);
xlabel('$\Lambda_{max} t$', 'interpreter', 'latex');
ylabel('$j$', 'interpreter', 'latex');
title('$\mathrm{Prediction}$', 'interpreter', 'latex', 'Fontsize', 14);
set(gca, 'Fontsize', 12)


subplot(3, 1, 3);
show = X_test(:, 1:t_contour) - X_pred(:, 1:t_contour);
[CR, h] = contourf(t(1:t_contour), x, show, 'linestyle', 'none');
hold on
plot(t_blow*ones(100, 1), linspace(1, n, 100), ':', 'linewidth', 2, 'color', 'black');
hold on
colormap(b2r(-20, 20));
cb3 = colorbar;
caxis([-20 20]);
pbaspect([4, 1, 1]);
xlim([0 16]);
ylim([1 n]);
xlabel('$\Lambda_{max} t$', 'interpreter', 'latex');
ylabel('$j$', 'interpreter', 'latex');
title('$\mathrm{Difference}$', 'interpreter', 'latex', 'Fontsize', 14);
set(gca, 'Fontsize', 12)






