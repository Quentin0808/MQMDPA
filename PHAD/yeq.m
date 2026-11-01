function [y, dy] = yeq(Par, npa)

    comps = Par.phases.Liquid.comps;
    Nc  = numel(comps);
    Npa = numel(npa);

    npa = npa(:);
    npa = max(npa, eps);
    nptot = sum(npa);

    if nptot <= 0
        error('yeq:InvalidNpa', 'Total number of pairs must be positive.');
    end

    y = zeros(Nc, 1);

    comp1 = Par.topo.comp1;
    comp2 = Par.topo.comp2;

    for j = 1:Npa
        i1 = comp1(j);
        i2 = comp2(j);

        if i1 == i2
            y(i1) = y(i1) + npa(j);
        else
            y(i1) = y(i1) + 0.5 * npa(j);
            y(i2) = y(i2) + 0.5 * npa(j);
        end
    end

    y = y / nptot;
    y = max(y, eps);

    if nargout > 1
        dy = zeros(Nc, Npa);

        for j = 1:Npa
            i1 = comp1(j);
            i2 = comp2(j);

            if i1 == i2
                dy(:, j) = -y / nptot;
                dy(i1, j) = dy(i1, j) + 1 / nptot;
            else
                dy(:, j) = -y / nptot;
                dy(i1, j) = dy(i1, j) + 0.5 / nptot;
                dy(i2, j) = dy(i2, j) + 0.5 / nptot;
            end
        end
    end
end