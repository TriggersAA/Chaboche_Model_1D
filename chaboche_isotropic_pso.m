
clear; clc;

%% --- Load Experimental Data ---
opts = detectImportOptions('1.5percent.csv', 'VariableNamingRule', 'preserve');
data = readtable('1.5percent.csv', opts);
N = 10099;  % Use N+1 points (2 to 1200)

total_strain = data{2:N+1, 1};    % Total strain
plastic_strain = data{2:N+1, 2};  % Plastic strain
stress_exp = data{2:N+1, 3};      % Experimental stress

total_strain2 = data{1:35, 1};  % 1st column
plastic_strain2 = data{1:35, 2};  % 2nd column
stress2 = data{1:35, 3};  % 3rd column

elastic_idx = total_strain2 >= 0;
strain_elastic = total_strain2(elastic_idx);
stress_elastic = stress2(elastic_idx);

% Estimate Young's modulus E
p = polyfit(total_strain2, stress2, 1);
E = p(1);

if length(strain_elastic) < 2
    error('Not enough elastic data points for E estimation');
end



%% --- Yield Stress (Fixed) ---
sigma_y0 = 1200; % MPa

% --- Chaboche model with isotropic hardening and beta-switching on X3 only ---
J2 = @(x) 1.5 * x.^2;

function error = chaboche_iso_hard(params, total_strain, plastic_strain, stress_exp, E, sigma_y0)
    J2 = @(x) 1.5 * x.^2;

    C1 = params(1); g1 = params(2);
    C2 = params(3); g2 = params(4);
    C3 = params(5); g31 = params(6); g32 = params(7);
    Q = params(8); b = params(9);

    N = length(total_strain);
    X1 = 0; X2 = 0; X3 = 0; R = 0;
    stress_model = zeros(N,1);

    for i = 2:N
        d_eps_p = plastic_strain(i) - plastic_strain(i-1);
        eps_e = total_strain(i) - plastic_strain(i);
        sigma_trial = E * eps_e;
        X = X1 + X2 + X3;
        s = sigma_trial - X;
        f = J2(s) - (sigma_y0 + R);

        if f > 0
            dp = d_eps_p;
            n_dir = sign(s);

            X1 = X1 + (2/3)*C1*dp*n_dir - g1 * X1 * dp;
            X2 = X2 + (2/3)*C2*dp*n_dir - g2 * X2 * dp;

            j2_x3 = J2(X3);
            x3_eq = abs(X3);
            if j2_x3 <= x3_eq
                beta = g31;
            else
                beta = g32 * (1 - x3_eq / j2_x3);
            end
            X3 = X3 + (2/3)*C3*dp*n_dir - beta * X3 * dp;

            R = R + b * (Q - R) * dp;
        end
        stress_model(i) = sigma_trial;
    end
    error = norm(stress_model - stress_exp);
end

% --- Optimization using Particle Swarm ---
x0 = [2000, 20, 1500, 15, 1000, 5, 80, 100, 10]; % Initial guess
lb = [1, 1, 1, 1, 1, 1, 1, 1, 0.001];
ub = [1e5, 1e3, 1e5, 1e3, 1e5, 1e2, 1e2, 1e3, 10];

obj_fun = @(x) chaboche_iso_hard(x, total_strain, plastic_strain, stress_exp, E, sigma_y0);

options = optimoptions('particleswarm','Display','iter','SwarmSize',100,'MaxIterations',200);
[x_opt, fval] = particleswarm(obj_fun, 9, lb, ub, options);

% --- Re-simulate with optimized parameters ---
C1 = x_opt(1); g1 = x_opt(2);
C2 = x_opt(3); g2 = x_opt(4);
C3 = x_opt(5); g31 = x_opt(6); g32 = x_opt(7);
Q = x_opt(8); b = x_opt(9);

X1 = 0; X2 = 0; X3 = 0; R = 0;
stress_model = zeros(N,1);

for i = 2:N
    d_eps_p = plastic_strain(i) - plastic_strain(i-1);
    eps_e = total_strain(i) - plastic_strain(i);
    sigma_trial = E * eps_e;
    X = X1 + X2 + X3;
    s = sigma_trial - X;
    f = J2(s) - (sigma_y0 + R);

    if f > 0
        dp = d_eps_p;
        n_dir = sign(s);

        X1 = X1 + (2/3)*C1*dp*n_dir - g1 * X1 * dp;
        X2 = X2 + (2/3)*C2*dp*n_dir - g2 * X2 * dp;

        j2_x3 = J2(X3);
        x3_eq = abs(X3);
        if j2_x3 <= x3_eq
            beta = g31;
        else
            beta = g32 * (1 - x3_eq / j2_x3);
        end
        X3 = X3 + (2/3)*C3*dp*n_dir - beta * X3 * dp;

        R = R + b * (Q - R) * dp;
    end
    stress_model(i) = sigma_trial;
end

% --- Plot Results ---
figure;
plot(total_strain, stress_exp, 'b', 'DisplayName', 'Experimental');
hold on;
plot(total_strain, stress_model, 'r--', 'DisplayName', 'Model Prediction');
xlabel('Total Strain'); ylabel('Stress (MPa)');
legend; grid on;
title('Chaboche Model with Isotropic Hardening (Particle Swarm)');

% --- Display optimized parameters ---
fprintf('Optimized Parameters:\n');
fprintf('C1  = %.4f\n', x_opt(1));
fprintf('g1  = %.4f\n', x_opt(2));
fprintf('C2  = %.4f\n', x_opt(3));
fprintf('g2  = %.4f\n', x_opt(4));
fprintf('C3  = %.4f\n', x_opt(5));
fprintf('g31 = %.4f\n', x_opt(6));
fprintf('g32 = %.4f\n', x_opt(7));
fprintf('Q   = %.4f\n', x_opt(8));  % Isotropic hardening saturation
fprintf('b   = %.4f\n', x_opt(9));  % Isotropic hardening rate