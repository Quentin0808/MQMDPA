function Par = build_comp_npa_matrix(Par)

    comps = Par.phases.Liquid.comps;
    ncomp = numel(comps);
    npatot = Par.topo.total_npa;

    A = zeros(ncomp, npatot);

    for j = 1:npatot
        key_base = Par.topo.base_key{j};
        localIdx = Par.topo.local_id_of_pair(j);

        im = Par.topo.comp1(j);
        in = Par.topo.comp2(j);

        Z = get_Z_by_localIdx(Par, key_base, localIdx);

        if ~isnumeric(Z) || numel(Z) ~= 2
            error('build_comp_npa_matrix:InvalidZ', ...
                'Invalid Z for pair "%s" (localIdx = %d).', ...
                key_base, localIdx);
        end

        Zm = Z(1);
        Zn = Z(2);

        if Zm == 0 || Zn == 0
            error('build_comp_npa_matrix:ZeroCoordinationNumber', ...
                'Zero coordination number found for pair "%s" (localIdx = %d).', ...
                key_base, localIdx);
        end

        A(im, j) = A(im, j) + 1 / Zm;
        A(in, j) = A(in, j) + 1 / Zn;
    end

    Par.idx.comp_npa_matrix = A;

end