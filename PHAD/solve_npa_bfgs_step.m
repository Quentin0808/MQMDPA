function [npa_new, lambda, B_new, info] = solve_npa_bfgs_step(Par, T, npa0, b, B, opts)

    if nargin < 6
        opts = struct();
    end

    if ~isfield(opts,'tol'); opts.tol = 1e-12; end
    if ~isfield(opts,'reg'); opts.reg = 1e-8; end
    if ~isfield(opts,'max_step_fraction'); opts.max_step_fraction = 0.8; end

    Aeq = Par.idx.comp_npa_matrix;
    b = b(:); % the system composition with Nc * 1 matrix 
    npa0 = npa0(:);
    Npa = length(npa0);

    npa0 = max(npa0, opts.tol);

    dg0 = grad_mqmdpa(Par, T, npa0);
    dg0 = dg0(:);

    % ============= B =============
    
    need_rebuild = isempty(B);

    if ~need_rebuild
        Bs = (B+B')/2;
        eigB = eig(Bs);
        need_rebuild = ...
        rcond(Bs) < 1e-30 || ...
        min(eigB) <= 0 || ...
        any(isnan(Bs(:))) || ...
        any(isinf(Bs(:)));
    end
     
    if need_rebuild
        B = Hessian(@(x) grad_mqmdpa(Par,T,x), npa0, 1e-6);
    end

    B = regularize_spd(B);
    
    % ------------------ B ------------------
    Bsolve = B + opts.reg * eye(Npa);
%       Bsolve = B;

    % ---------------- cG, cA ---------------
    [cG, cU, info_ci] = solve_ci_coefficients(Bsolve, dg0, Aeq.');

    % ------------- Calculate lambda ------------
    K   = Aeq * cU;
    rhs = b - Aeq*npa0 - Aeq*cG;
    lambda = K \ rhs;
    dnpa_full = cG + cU * lambda;

    % npa >= tol
    idx = dnpa_full < 0;
    if any(idx)
        alpha_max = min((npa0(idx) - opts.tol) ./ (-dnpa_full(idx)));
        alpha = min(1.0, opts.max_step_fraction * alpha_max);
    else
        alpha = 1;
    end

    if alpha < 0 || isnan(alpha) || isinf(alpha)
        alpha = 0;
    end

    dnpa = alpha * dnpa_full;
    npa_new = npa0 + dnpa;
    npa_new = max(npa_new, opts.tol);

    % ===== BFGS =====
    dg1 = grad_mqmdpa(Par, T, npa_new);
    dg1 = dg1(:);

    s = npa_new - npa0;
    y = dg1 - dg0;

    ys = y.' * s;

    if ys > 1e-12 * norm(y) * norm(s)
        Bs = B * s;
        sBs = s.' * Bs;

        if sBs > 0
            B_new = B + (y*y.')/ys - (Bs*Bs.')/sBs;
        else
            B_new = B;
        end
    else
        B_new = B;
    end

    % =========== symmetrize ===========
    
    B_new = 0.5 * (B_new + B_new.');

    % =========== diagnostics ==========
    info = struct();
    info.alpha = alpha;
    info.dg0 = dg0;
    info.dg1 = dg1;
    info.lambda = lambda;
    info.cG = cG;
    info.cU = cU;
    info.info_ci = info_ci;
    info.dnpa_full = dnpa_full;
    info.dnpa = dnpa;
    info.step_norm = norm(dnpa);
    info.mass_error_before = norm(Aeq*npa0 - b);
    info.mass_error_after  = norm(Aeq*npa_new - b);
    info.min_npa_after = min(npa_new);
    info.kkt_error = norm(dg1 - Aeq.'*lambda);
    info.ys = ys;
    info.cond_B = cond(Bsolve);
    info.cond_K = cond(K);
end