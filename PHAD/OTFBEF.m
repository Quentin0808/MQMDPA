function [gb, dgb] = OTFBEF(Par, mqstate, idp)

    tol = 1e-12;
    
    T = mqstate.T;
    npa = real(mqstate.npa(:));
    y = real(mqstate.y(:));
    dy = real(mqstate.dy);
    
    npa = max(npa, tol);
    y = max(y, tol);
    
    Npa   = Par.topo.total_npa;
    phase = Par.phases.Liquid;
    comps = phase.comps;

    if Npa ~= numel(npa)
        error('OTFBEF:SizeMismatch', ...
            'Number of pairs mismatch: expected %d, got %d.', ...
            Npa, numel(npa));
    end

    gb  = 0;
    dgb = zeros(Npa,1);

    if Par.topo.is_homo(idp)
        return;
    end

    localIdx = Par.topo.local_id_of_pair(idp);

    mi = Par.topo.comp1(idp);
    ni = Par.topo.comp2(idp);

    m = comps{mi};
    n = comps{ni};

    keymn = canonical_key(Par, m, n);

    bmn = Par.topo.base_id_of_pair(idp);
    ids_mn = Par.idx.npa_ids_of_base{bmn};

    id_mm = Par.idx.homo_ids(mi);
    id_nn = Par.idx.homo_ids(ni);

    bmm = Par.topo.base_id_of_pair(id_mm);
    bnn = Par.topo.base_id_of_pair(id_nn);

    ids_mm = Par.idx.npa_ids_of_base{bmm};
    ids_nn = Par.idx.npa_ids_of_base{bnn};

    nmm = sum(npa(ids_mm));
    nnn = sum(npa(ids_nn));
    nmn = sum(npa(ids_mn));

    yp  = 0;
    py  = 0;
    dyp = zeros(Npa,1);
    dpy = zeros(Npa,1);

    gb_bin  = 0;
    gb_ten  = 0;
    dgb_bin = zeros(Npa,1);
    dgb_ten = zeros(Npa,1);

    fmn = nmm;
    fnm = nnn;
    wmn = nmm + nnn + nmn;

    dfmn = zeros(Npa,1);
    dfnm = zeros(Npa,1);
    dwmn = zeros(Npa,1);

    dfmn(ids_mm) = 1;
    dfnm(ids_nn) = 1;

    dwmn(ids_mm) = 1;
    dwmn(ids_nn) = 1;
    dwmn(ids_mn) = 1;

    if isfield(Par.topo, 'others_of_pair')
        others = Par.topo.others_of_pair{idp};
    else
        error('OTFBEF:MissingTopologyCache', ...
            'Please call build_others_of_pair(Par) before OTFBEF.');
    end

    %--------------------------------------------------------------
    % Build fmn, fnm, wmn from preprocessed other components
    %--------------------------------------------------------------
    for a = 1:numel(others)

        si   = others(a).si;
        Dmns = others(a).D;

        ids_ss = others(a).ids_ss;
        ids_ms = others(a).ids_ms;
        ids_ns = others(a).ids_ns;

        nss = sum(npa(ids_ss));
        nms = sum(npa(ids_ms));
        nns = sum(npa(ids_ns));

        D1 = Dmns(1);
        D2 = Dmns(2);
        DS = D1 + D2;

        yp  = yp  + y(si) * D1;
        py  = py  + y(si) * D2;

        dyp = dyp + dy(si,:).' * D1;
        dpy = dpy + dy(si,:).' * D2;

        fmn = fmn + nss * D1^2 + D1 * nms;
        fnm = fnm + nss * D2^2 + D2 * nns;
        wmn = wmn + nss * DS^2 + DS * (nms + nns);

        dfmn(ids_ss) = D1^2;
        dfnm(ids_ss) = D2^2;
        dwmn(ids_ss) = DS^2;

        dfmn(ids_ms) = D1;
        dfnm(ids_ns) = D2;

        dwmn(ids_ms) = DS;
        dwmn(ids_ns) = DS;

        for b = 1:numel(others(a).cross)

            Dmnt   = others(a).cross(b).D;
            ids_st = others(a).cross(b).ids_st;

            nst = sum(npa(ids_st));

            E1 = Dmnt(1);
            E2 = Dmnt(2);
            ES = E1 + E2;

            fmn = fmn + D1 * E1 * nst;
            fnm = fnm + D2 * E2 * nst;
            wmn = wmn + DS * ES * nst;

            dfmn(ids_st) = D1 * E1;
            dfnm(ids_st) = D2 * E2;
            dwmn(ids_st) = DS * ES;
        end
    end

    if ~isreal(wmn) || abs(wmn) < tol
        return;
    end

    kapmn = real(fmn / wmn);
    kapnm = real(fnm / wmn);
    kapmn = max(kapmn, tol);
    kapnm = max(kapnm, tol);

    dkapmn = (wmn * dfmn - fmn * dwmn) / wmn^2;
    dkapnm = (wmn * dfnm - fnm * dwmn) / wmn^2;

    e_idp = zeros(Npa,1);
    e_idp(idp) = 1;

    if npa(idp) > eps
        dlog_npa = e_idp / npa(idp);
    else
        dlog_npa = zeros(Npa,1);
    end

    if abs(kapmn) > eps
        dlog_kapmn = dkapmn / kapmn;
    else
        dlog_kapmn = zeros(Npa,1);
    end

    if abs(kapnm) > eps
        dlog_kapnm = dkapnm / kapnm;
    else
        dlog_kapnm = zeros(Npa,1);
    end

    %--------------------------------------------------------------
    % Binary contribution
    %--------------------------------------------------------------
    if isfield(phase, keymn)

        entry = phase.(keymn);

        if localIdx <= numel(entry) && isfield(entry(localIdx), 'dg_list')

            lst = entry(localIdx).dg_list;

            for r = 1:size(lst,1)

                i = lst(r,1);
                j = lst(r,2);

                gij = lst(r,3) + lst(r,4) * T;

                gmnij = gij * kapmn^i * kapnm^j;

                gb_bin = gb_bin + gmnij;

                dgb_bin = dgb_bin + gmnij * ...
                    (dlog_npa + i*dlog_kapmn + j*dlog_kapnm);
            end
        end
    end

    %--------------------------------------------------------------
    % Ternary contribution
    %--------------------------------------------------------------
    epmn  = real(y(mi) + yp);
    epnm  = real(y(ni) + py);

    depmn = dy(mi,:).' + dyp;
    depnm = dy(ni,:).' + dpy;

    if epmn > tol
        epmn_safe  = epmn;
        depmn_safe = depmn;
    else
        epmn_safe  = tol;
        depmn_safe = zeros(Npa,1);
    end

    if epnm > tol
        epnm_safe  = epnm;
        depnm_safe = depnm;
    else
        epnm_safe  = tol;
        depnm_safe = zeros(Npa,1);
    end

    rest   = max(real(1 - epmn_safe - epnm_safe), tol);
    termKm = max(real(1 - y(mi) / epmn_safe), tol);
    termKn = max(real(1 - y(ni) / epnm_safe), tol);

    for a = 1:numel(others)

        si    = others(a).si;
        Dmns  = others(a).D;
        keymns = others(a).keymns;

        if ~isfield(phase, keymns)
            continue;
        end

        entry = phase.(keymns);

        if localIdx > numel(entry)
            continue;
        end

        item = entry(localIdx);

        if ~isfield(item, 'dg_list') || isempty(item.dg_list)
            continue;
        end

        lst = item.dg_list;

        D1 = Dmns(1);
        D2 = Dmns(2);

        for r = 1:size(lst,1)

            i = lst(r,1);
            j = lst(r,2);
            k = lst(r,3);
            aa = lst(r,4);
            bb = lst(r,5);

            if k < 1
                continue;
            end

            gijk = aa + bb*T;

            toop  = (1 - D1 - D2) * rest^(k-1);
            kohlm = (D1 / epmn_safe) * termKm^(k-1);
            kohln = (D2 / epnm_safe) * termKn^(k-1);

            denom = toop + kohlm + kohln;

            denom = real(denom);

            if abs(denom) < tol
                continue;
            end         
            
            if y(si) > tol
                dlog_y = dy(si,:).' / y(si);
            else
                dlog_y = zeros(Npa,1);
            end

            if k == 1
                dtoop = zeros(Npa,1);
            else
                dtoop = -(1 - D1 - D2) * (k-1) * rest^(k-2) ...
                      * (depmn_safe + depnm_safe);
            end

            if k == 1
                dkohlm = -D1 / epmn_safe^2 * depmn_safe;
                dkohln = -D2 / epnm_safe^2 * depnm_safe;
            else
                dkohlm = -D1 / epmn_safe^2 * termKm^(k-2) ...
                       * (termKm * depmn_safe ...
                       + (k-1) * (dy(mi,:).' * epmn_safe ...
                       - y(mi) * depmn_safe) / epmn_safe);

                dkohln = -D2 / epnm_safe^2 * termKn^(k-2) ...
                       * (termKn * depnm_safe ...
                       + (k-1) * (dy(ni,:).' * epnm_safe ...
                       - y(ni) * depnm_safe) / epnm_safe);
            end

            gmnsijk = gijk * kapmn^i * kapnm^j * y(si) * denom;

            gb_ten = gb_ten + gmnsijk;

            dgb_ten = dgb_ten + gmnsijk * ...
                (dlog_npa + i*dlog_kapmn + j*dlog_kapnm ...
                + dlog_y + (dtoop + dkohlm + dkohln) / denom);
        end
    end

    gb  = real(gb_bin + gb_ten);
    dgb = real(dgb_bin + dgb_ten);
    dgb = dgb(:);

    if ~isfinite(gb) || any(~isfinite(dgb))
        error('OTFBEF:NonFiniteOutput', 'gb or dgb contains NaN/Inf.');
    end
end