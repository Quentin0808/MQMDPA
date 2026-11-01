function Par = build_others_of_pair(Par)

    phase = Par.phases.Liquid;
    comps = phase.comps;
    Nc    = Par.Nc;
    Npa   = Par.topo.total_npa;

    Par.topo.others_of_pair = cell(Npa,1); % centers on the selected pair which has Npa ways of choice. 

    for idp = 1:Npa

        if Par.topo.is_homo(idp)
            Par.topo.others_of_pair{idp} = [];
            continue;
        end

        mi = Par.topo.comp1(idp);
        ni = Par.topo.comp2(idp);

        m = comps{mi};
        n = comps{ni};

        keymn = canonical_key(Par, m, n);

        other_ids = setdiff(1:Nc, [mi ni], 'stable');
        
        others = struct([]);

        for a = 1:numel(other_ids)

            si = other_ids(a);
            s  = comps{si};

            keyss  = [s s];
            keyms  = canonical_key(Par, m, s);
            keyns  = canonical_key(Par, n, s);
            keymns = [keymn s];

            Dmns = get_nested_field(phase, 'D', keymns, [0 0]);

            bss = Par.idx.base_key_to_id.(keyss);
            bms = Par.idx.base_key_to_id.(keyms);
            bns = Par.idx.base_key_to_id.(keyns);

            others(a).si = si;
            others(a).D  = Dmns;

            others(a).ids_ss = Par.idx.npa_ids_of_base{bss};
            others(a).ids_ms = Par.idx.npa_ids_of_base{bms};
            others(a).ids_ns = Par.idx.npa_ids_of_base{bns};

            others(a).keymns = keymns;
            others(a).cross  = struct([]);

            for b = a+1:numel(other_ids)

                ti = other_ids(b);
                t  = comps{ti};

                keyst  = canonical_key(Par, s, t);
                keymnt = [keymn t];

                Dmnt = get_nested_field(phase, 'D', keymnt, [0 0]);

                bst = Par.idx.base_key_to_id.(keyst);

                cidx = b - a;

                others(a).cross(cidx).ti = ti;
                others(a).cross(cidx).D  = Dmnt;
                others(a).cross(cidx).ids_st = Par.idx.npa_ids_of_base{bst};
            end
        end

        Par.topo.others_of_pair{idp} = others;
    end
end