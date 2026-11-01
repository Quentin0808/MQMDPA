function key = canonical_key(Par, varargin)

    comps = Par.phases.Liquid.comps;
    n = numel(varargin);

    if n < 2 || n > 3
        error('canonical_key:InvalidInput', ...
              'Function requires 2 or 3 components.');
    end

    idx = zeros(1,n);

    % ?????? comps ????
    for i = 1:n
        idx(i) = find(strcmp(comps,varargin{i}),1);

        if isempty(idx(i))
            error('canonical_key:UnknownComponent', ...
                'Component "%s" not found.', varargin{i});
        end
    end

    % ? comps ????
    [~,order] = sort(idx);

    sorted = varargin(order);

    % ?????
   key = [sorted{:} ];

end