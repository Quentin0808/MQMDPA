function state = build_state(Par, T, npa)

    npa = npa(:);

    Npa = Par.topo.total_npa;
    Nc = Par.Nc;

    if numel(npa) ~= Npa
        error('build_state:SizeMismatch', ...
            'Expected npa size %d, got %d.', Npa, numel(npa));
    end

    state.T = T;
    state.npa = npa;

    % component amounts
    A = Par.idx.comp_npa_matrix;
    state.ncomp = A * npa;

    state.nctot = sum(state.ncomp);
    state.nptot = sum(npa);

    eps0 = 1e-15;

    if state.nctot <= 0
        error('build_state:InvalidState','Total component amount <= 0.');
    end

    if state.nptot <= 0
        error('build_state:InvalidState','Total pair amount <= 0.');
    end

    state.xc  = state.ncomp / state.nctot;
    state.xpa = npa / state.nptot;

    % safe versions
    state.xc_safe  = max(state.xc,  eps0);
    state.xpa_safe = max(state.xpa, eps0);

    % equivalent coordination fractions
    [state.y, state.dy] = yeq(Par, npa);
    state.y_safe = max(state.y, eps0);

    % pure G0
    comps = Par.phases.Liquid.comps;
    state.G0 = zeros(Nc,1);

    for i = 1:Nc
        state.G0(i) = G0_pure(Par, comps{i}, 'LIQ', T);
    end
end