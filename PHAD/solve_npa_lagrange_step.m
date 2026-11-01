function [npa_new, mu, info] = solve_npa_lagrange_step(Par, T, npa0, b, opts)
% Solve one constrained Newton/Lagrange correction step:
%
%   minimize G(npa)
%   subject to:
%       Aeq * npa = b
%       npa >= tol
%
% Uses:
%   dnpa = cG + cA * mu
%
% where:
%   H*cG = -gradG
%   H*cA = Aeq'
%
% Then enforce:
%   Aeq*(npa0 + dnpa) = b

    if nargin < 5
        opts = struct();
    end

    if ~isfield(opts,'h')
        opts.h = 1e-5;
    end

    if ~isfield(opts,'tol')
        opts.tol = 1e-12;
    end

    if ~isfield(opts,'reg')
        opts.reg = 1e-8;
    end

    if ~isfield(opts,'max_step_fraction')
        opts.max_step_fraction = 0.8;
    end

    if ~isfield(opts,'use_regularization')
        opts.use_regularization = true;
    end

    Aeq = Par.idx.comp_npa_matrix;
    b = b(:);
    npa0 = npa0(:);

    Npa = length(npa0);

    % Safety: keep npa0 positive
    npa0 = max(npa0, opts.tol);

    % Gradient and Hessian
    dg = grad_mqmdpa(Par, T, npa0);
    dg = dg(:);

    H = Hessian(@(x) grad_mqmdpa(Par, T, x), npa0, opts.h);

    % Symmetrize Hessian
    if norm(H - H.', 'fro') > 1e-10
        warning('H is not symmetric');
    end

    % Regularization for ill-conditioned Hessian
    if opts.use_regularization
        Hsolve = H + opts.reg * eye(Npa);
    else
        Hsolve = H;
    end

    % Phase matrix solve:
    %
    % H*cG = -dg
    % H*cA = Aeq'
    [cG, cA, info_ci] = solve_ci_coefficients(Hsolve, dg, Aeq.');
    
    
    % Equilibrium matrix:
    %
    % Aeq*cA*mu = b - Aeq*npa0 - Aeq*cG
    K = Aeq * cA;
    rhs = b - Aeq*npa0 - Aeq*cG;
    mu = K \ rhs;
    
    % Full Newton correction
    dnpa_full = cG + cA * mu;

    % Step damping to keep npa >= tol
    idx = dnpa_full < 0;

    if any(idx)
        alpha_max = min((npa0(idx) - opts.tol) ./ (-dnpa_full(idx)));
        alpha = min(1.0, opts.max_step_fraction * alpha_max);
    else
        alpha = 1.0;
    end

    if alpha < 0 || isnan(alpha) || isinf(alpha)
        alpha = 0;
    end

    dnpa = alpha * dnpa_full;
    npa_new = npa0 + dnpa;

    % Final safety
    npa_new = max(npa_new, opts.tol);

    % Diagnostics
    info = struct();
    info.H = H;
    info.Hsolve = Hsolve;
    info.dg = dg;
    info.cG = cG;
    info.cA = cA;
    info.K = K;
    info.rhs = rhs;
    info.dnpa_full = dnpa_full;
    info.dnpa = dnpa;
    info.alpha = alpha;
    info.mass_error_before = norm(Aeq*npa0 - b);
    info.mass_error_after = norm(Aeq*npa_new - b);
    info.min_npa_before = min(npa0);
    info.min_npa_after = min(npa_new);
    info.cond_H = cond(Hsolve);
    info.cond_K = cond(K);
    info.step_norm = norm(dnpa);
    info.kkt_error = norm(dg - Aeq.'*mu);
    info.info_ci = info_ci;
end