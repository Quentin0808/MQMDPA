clc; clear; close all;

format long

% % ============================================================
% % Build parameters
% % ============================================================

Par = params_NaMgCl();

Par = build_pair_topology(Par);
Par = build_comp_npa_matrix(Par);
Par = build_others_of_pair(Par);

Aeq = Par.idx.comp_npa_matrix;
nc    = Par.Nc;
nptot = Par.topo.total_npa;

T = 1273;

% ============================================================
% Options
% ============================================================

opts = struct();
opts.tol = 1e-12;
opts.reg = 1e-8;
opts.max_step_fraction = 0.8;

% ============================================================
% Composition scan
% ============================================================

NUM = 100;

xlist = linspace(1e-20, 1.0, NUM);

Glist = zeros(NUM,1);
npa_store = zeros(nptot, NUM);
lambda_store  = zeros(nc, NUM);

mass_error = zeros(NUM,1);
step_error = zeros(NUM,1);
kkt_error  = zeros(NUM,1);
iter_store = zeros(NUM,1);

% ============================================================
% Initial composition
% ============================================================

% intial npa for pair numbers

  npa0 = rand(nptot,1);
  npa = npa0/sum(npa0)*3;

% initial B

  B = eye(nptot);
  
% ============================================================
% Main loop
% ============================================================

tic

for i = 1:NUM

    x = xlist(i);
%     b = [x; 0.5*(1-x); 0.25*(1-x); 0.25*(1-x)];
    b = [1-x; x];
% make sure current npa is real and positive
    npa = real(npa);
    npa = max(npa, opts.tol);
%     
    if i <= 5 || i >= 95
        maxiter = 200;
    else
        maxiter = 30;
    end

    for iter = 1:maxiter
                
        [npa_new, lambda, B, info] = solve_npa_bfgs_step(Par, T, npa, b, B, opts);
        
        % safety check
        if any(~isreal(npa_new))
            error('Complex npa detected at i = %d, iter = %d', i, iter);
        end

        if any(~isfinite(npa_new))
            error('NaN/Inf npa detected at i = %d, iter = %d', i, iter);
        end
                    
        npa = real(npa_new);
        npa = max(npa, opts.tol);

        if info.mass_error_after < 1e-8 && ...
           info.step_norm < 1e-7 && ...
           info.kkt_error < 1e-3
            break
        end
    end

    Glist(i) = G_mqmdpa(Par, T, npa);

    npa_store(:,i) = npa;
    lambda_store(:,i)  = lambda;

    mass_error(i) = info.mass_error_after;
    step_error(i) = info.step_norm;
    kkt_error(i)  = info.kkt_error;
    iter_store(i) = iter;

    fprintf(['i = %d/%d, x = %.5f, G = %.8e, iter = %d, ', ...
             'mass = %.3e, step = %.3e, kkt = %.3e\n'], ...
             i, NUM, x, Glist(i), iter, ...
             mass_error(i), step_error(i), kkt_error(i));
end

toc

% ============================================================
% Plot Gibbs energy  
% ============================================================

figure;
plot(xlist, Glist, 'LineWidth', 2)
xlabel('Composition');
ylabel('Gibbs Energy');
title('G vs Composition by BFGS solver');
grid on;


% ============================================================
% Plot diagnostics
% ============================================================
figure;
semilogy(xlist, mass_error, 'o-', 'LineWidth', 1.5); hold on;
semilogy(xlist, step_error, 's-', 'LineWidth', 1.5);
semilogy(xlist, kkt_error, '^-', 'LineWidth', 1.5);
legend('mass error','step norm','KKT error');
xlabel('Composition');
ylabel('Error');
title('Solver diagnostics');
grid on;

% ============================================================
% Plot pair amounts
% ============================================================
figure;
plot(xlist, npa_store.', 'LineWidth', 1.5);
xlabel('Composition');
ylabel('npa');
title('Pair populations');
grid on;

% ============================================================
% Plot chemical potential
% ============================================================
figure;
plot(xlist, lambda_store.', 'LineWidth', 1.5);
xlabel('Composition');
ylabel('lambda');
title('Chemical potential');
grid on;