function [A, B, phi, Lambda, U_tau, U_1, U_2] = HDMDC_Execution(data, control, p, r, q, tau)

n_data = size(data, 1);
m      = size(data, 2);

[X_delayed, Y_delayed] = Hankel_DMD_Matrices(data, n_data, m, q);

n_control = size(control, 1);

control_delayed = Hankel_DMD_Matrices(control, n_control, m, q);

C = [X_delayed; control_delayed];

[U, S, V] = svd(C, 'econ');

U = U(1:end, 1:p);
S = S(1:p, 1:p);
V = V(1:end, 1:p);

U_1 = U(1:n_data*q, :);
U_2 = U(n_data*q+1:end, :);

[U_tau, S_tau, V_tau] = svd(Y_delayed, 'econ');

U_tau = U_tau(1:end, 1:r);
% S_tau = S_tau(1:r, 1:r);
% V_tau = V_tau(1:end, 1:r);

A = U_tau'*Y_delayed*(V/S)*U_1'*U_tau;
B = U_tau'*Y_delayed*(V/S)*U_2';

[W, Lambda] = eig(A);
phi = Y_delayed*(V/S)*U_1'*U_tau*W;

[Lambda, I] = sort(abs(diag(Lambda)), 'descend');
phi = phi(:, I);

phi_n = zeros(size(phi, 1)/q, size(phi, 2));
for(i = 1:size(phi, 2))
    phi_n(:, i) = phi(end-size(phi, 1)/q+1:end, i);
end
phi = phi_n;

return 
