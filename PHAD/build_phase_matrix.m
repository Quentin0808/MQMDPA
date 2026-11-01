function phase = build_phase_matrix(phase, T)
    
    y = phase.var;                 % ???????MQMDPA ?? npa?CEF ?? y
    A = phase.comp_matrix;          % M = A*y

    dg = phase.grad(T, y);

    if isfield(phase, 'hess') && ~isempty(phase.hess)
        H = phase.hess(T, y);
    else
        H = Hessian(@(x) phase.grad(T,x), y, phase.hess_step);
    end

    H = regularize_spd(H);

    [cG, cA, info_ci] = solve_ci_coefficients(H, dg, A.');

    phase.cG = cG;
    phase.cA = cA;
    phase.M  = A * y;
    phase.Gm = phase.eval(T, y);

    phase.dM_dmu = A * cA;
    phase.info_ci = info_ci;
end