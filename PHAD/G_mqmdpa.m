function [G, dG]= G_mqmdpa(Par, T, npa)
 
    state = build_state(Par, T, npa);

    R = Par.const.R;
    RT = R * T;
    npa = state.npa;
    Nc = Par.Nc;
    Npa = Par.topo.total_npa;
    ncomp = state.ncomp;
    xc = state.xc;
    xpa = state.xpa;
    
    y = state.y;
  
    G = 0;
    G0 = state.G0;
    gradGexs = zeros(Npa, 1);
    
    %----------------------------------------------------------
    % Initialize contributions
    %----------------------------------------------------------
    
    eps0 = 1e-15;

    for i = 1:Nc
        
        G = G + ncomp(i) * G0(i);

        if xc(i) > eps0
            G = G + RT * ncomp(i) * log(xc(i));
        end
        
        id_homo = Par.idx.homo_ids(i);
        
        if npa(id_homo) > eps0 && y(i) > eps0
            G = G + RT * npa(id_homo) * log(xpa(id_homo) / y(i)^2);
        end
    end

    %----------------------------------------------------------
    % Hetero pair configurational entropy + excess term
    %----------------------------------------------------------
    for b = 1:Par.topo.nBasePairs
        ids = Par.idx.npa_ids_of_base{b};

        % skip homo base pairs
        if Par.topo.is_homo(ids(1))
            continue;
        end

        i = Par.topo.comp1(ids(1));
        j = Par.topo.comp2(ids(1));

        q = numel(ids);

        for k = ids
            if npa(k) > eps0 && y(i) > eps0 && y(j) > eps0
                G = G + RT * npa(k) * log(q * xpa(k) / (2 * y(i) * y(j)));
            end
            
           [gex, dgex] = OTFBEF(Par, state, k);
           
            G = G + npa(k) * gex;
            gradGexs = gradGexs + npa(k) * dgex;
        end
    end

    if nargout > 1

    % ----- topology arrays -----
        comp1   = Par.topo.comp1;
        comp2   = Par.topo.comp2;
        is_homo = Par.topo.is_homo;
        base_id = Par.topo.base_id_of_pair;

        Z1_all  = Par.topo.Z1;
        Z2_all  = Par.topo.Z2;

        qpair   = Par.topo.base_count(base_id);

    % ----- pair derivative without string lookup -----
        gpa = zeros(Npa, 1);

        for k = 1:Npa

            i = comp1(k);
            j = comp2(k);

            Z1 = Z1_all(k);
            Z2 = Z2_all(k);

            if is_homo(k)

                gpa(k) = 2/Z1 * G0(i) + RT * (2/Z1 * log(xc(i)) ...
                       + log(xpa(k) / y(i)^2));

            else

                gpa(k) = 1/Z1 * G0(i) + 1/Z2 * G0(j) ...
                       + RT * (1/Z1 * log(xc(i)) + 1/Z2 * log(xc(j)) ...
                       + log(qpair(k) * xpa(k) / (2 * y(i) * y(j))));
            end
        end
        
        dG = gpa + gradGexs;
        
    end
end