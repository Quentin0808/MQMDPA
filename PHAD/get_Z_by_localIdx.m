function Z = get_Z_by_localIdx(Par,key_base,localIdx)

    Z = [6 6];

    if nargin < 3 || isempty(localIdx)
        localIdx = 1;
    end

    if ~isfield(Par.phases.Liquid,'Z')
        return
    end

    Zstruct = Par.phases.Liquid.Z;

    if isfield(Zstruct,key_base)

        Zdata = Zstruct.(key_base);

        if isnumeric(Zdata) && numel(Zdata)==2
            Z = Zdata;
            return
        end

    end

    key_full = sprintf('%s%d',key_base,localIdx);

    if isfield(Zstruct,key_full)

        Zdata = Zstruct.(key_full);

        if isnumeric(Zdata) && numel(Zdata)==2
            Z = Zdata;
            return
        end

    end

end