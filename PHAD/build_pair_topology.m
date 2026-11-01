function Par = build_pair_topology(Par)

    comps = Par.phases.Liquid.comps;
    Nc = numel(comps);

    %--------------------------------------------------------------
    % Step 1. Count total distinguishable pairs (total_npa)
    %--------------------------------------------------------------
    total_npa = Nc;   % homo pairs: A-A, B-B, ...

    for i = 1:Nc
        Ai = comps{i};
        for j = i+1:Nc
            Bj = comps{j};
            key = [Ai Bj];

            if isfield(Par.phases.Liquid, 'qsro') && isfield(Par.phases.Liquid.qsro, key)
                Q = Par.phases.Liquid.qsro.(key);
            else
                Q = 1;
            end

            total_npa = total_npa + Q;
        end
    end

    % number of base pairs = homo + hetero combinations
    nBasePairs = Nc + Nc*(Nc-1)/2;

    %--------------------------------------------------------------
    % Step 2. Preallocate memory for pair-level arrays
    %--------------------------------------------------------------
    pair_name = cell(total_npa, 1);         % e.g. 'A-B#1'
    pair_key = cell(total_npa, 1);         % e.g. 'AB'
    base_key = cell(total_npa, 1);         % e.g. 'AB'

    comp1 = zeros(total_npa, 1);            % first component index
    comp2 = zeros(total_npa, 1);            % second component index

    is_homo = false(total_npa, 1);          % logical vector
    base_id_of_pair = zeros(total_npa, 1); % which base pair this pair belongs to
    local_id_of_pair = zeros(total_npa, 1); % local distinguishable index within same base pair

    %--------------------------------------------------------------
    % Step 3. Preallocate memory for base-pair-level arrays
    %--------------------------------------------------------------
    base_pair_name = cell(nBasePairs, 1);   % e.g. 'A-B'
    base_pair_key = cell(nBasePairs, 1);   % e.g. 'AB'

    base_first = zeros(nBasePairs, 1);      % first global pair id in this base pair
    base_last = zeros(nBasePairs, 1);      % last  global pair id in this base pair
    base_count = zeros(nBasePairs, 1);      % number of distinguishable pairs in this base pair

    % optional fast lookup map: key -> base_id
    base_key_to_id = struct();

    %--------------------------------------------------------------
    % Step 4. Fill homo pairs
    %--------------------------------------------------------------
    idp = 0;       % global pair id
    base_id = 0;   % base pair id

    for i = 1:Nc
        Ai  = comps{i};
        key = [Ai Ai];

        idp = idp + 1;
        base_id = base_id + 1;

        pair_name{idp} = sprintf('%s-%s', Ai, Ai);
        pair_key{idp}  = key;
        base_key{idp}  = key;

        comp1(idp) = i;
        comp2(idp) = i;

        is_homo(idp) = true;
        base_id_of_pair(idp) = base_id;
        local_id_of_pair(idp) = 1;

        base_pair_name{base_id} = sprintf('%s-%s', Ai, Ai);
        base_pair_key{base_id} = key;
        base_first(base_id) = idp;
        base_last(base_id) = idp;
        base_count(base_id) = 1;

        base_key_to_id.(key) = base_id;
    end

    %--------------------------------------------------------------
    % Step 5. Fill hetero distinguishable pairs
    %--------------------------------------------------------------
    for i = 1:Nc
        Ai = comps{i};

        for j = i+1:Nc
            Bj  = comps{j};
            key = [Ai Bj];

            if isfield(Par.phases.Liquid, 'qsro') && isfield(Par.phases.Liquid.qsro, key)
                Q = Par.phases.Liquid.qsro.(key);
            else
                Q = 1;
            end

            base_id = base_id + 1;

            base_pair_name{base_id} = sprintf('%s-%s', Ai, Bj);
            base_pair_key{base_id}  = key;
            base_key_to_id.(key) = base_id;

            first_id = idp + 1;

            for k = 1:Q
                idp = idp + 1;

                if Q == 1
                    pair_name{idp} = sprintf('%s-%s', Ai, Bj);
                else
                    pair_name{idp} = sprintf('%s-%s#%d', Ai, Bj, k);
                end

                pair_key{idp}  = key;
                base_key{idp}  = key;

                comp1(idp) = i;
                comp2(idp) = j;

                is_homo(idp) = false;
                base_id_of_pair(idp)  = base_id;
                local_id_of_pair(idp) = k;
            end

            last_id = idp;

            base_first(base_id) = first_id;
            base_last(base_id)  = last_id;
            base_count(base_id) = Q;
        end
    end

    %--------------------------------------------------------------
    % Step 6. Safety check
    %--------------------------------------------------------------
    if idp ~= total_npa
        error('build_pair_topology: total_npa mismatch.');
    end

    if base_id ~= nBasePairs
        error('build_pair_topology: nBasePairs mismatch.');
    end

    %--------------------------------------------------------------
    % Step 7. Build quick index lists
    %--------------------------------------------------------------
    npa_ids_of_base = cell(nBasePairs, 1);
    for b = 1:nBasePairs
        npa_ids_of_base{b} = base_first(b):base_last(b);
    end

    homo_ids   = find(is_homo);
    hetero_ids = find(~is_homo);
    
    %--------------------------------------------------------------
    % Step 8. Build coordination numbers of pairs
    %--------------------------------------------------------------    
    Z1_of_pair = zeros(total_npa,1);
    Z2_of_pair = zeros(total_npa,1);

    for k = 1:total_npa
        key = pair_key{k};
        lid = local_id_of_pair(k);
        Z = get_Z_by_localIdx(Par, key, lid);

        Z1_of_pair(k) = Z(1);
        Z2_of_pair(k) = Z(2);
    end

    %--------------------------------------------------------------
    % Step 9. Save to Par
    %--------------------------------------------------------------
    Par.Nc = Nc;

    % pair topology
    Par.topo.total_npa = total_npa;
    Par.topo.pair_name = pair_name;
    Par.topo.pair_key = pair_key;
    Par.topo.base_key = base_key;

    Par.topo.comp1 = comp1;
    Par.topo.comp2 = comp2;

    Par.topo.is_homo = is_homo;
    Par.topo.base_id_of_pair  = base_id_of_pair;
    Par.topo.local_id_of_pair = local_id_of_pair;

    % base pair topology
    Par.topo.nBasePairs = nBasePairs;
    Par.topo.base_pair_name = base_pair_name;
    Par.topo.base_pair_key = base_pair_key;
    Par.topo.base_first = base_first;
    Par.topo.base_last = base_last;
    Par.topo.base_count = base_count;
    Par.topo.Z1 = Z1_of_pair;
    Par.topo.Z2 = Z2_of_pair;

    % quick indices
    Par.idx.base_key_to_id = base_key_to_id;
    Par.idx.npa_ids_of_base = npa_ids_of_base;
    Par.idx.homo_ids = homo_ids;
    Par.idx.hetero_ids = hetero_ids;
end