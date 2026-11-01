function B = regularize_spd(B)

    % ????
    B = (B + B') / 2;

    % ?????
    [V,D] = eig(B);

    d = diag(D);

    % ?????
    dmax = max(abs(d));

    if isempty(dmax) || dmax == 0 || isnan(dmax) || isinf(dmax)
        dmax = 1;
    end

    % ???????
    floorEig = max(1e-10 * dmax, 1e-12);

    % clipping
    d(d < floorEig) = floorEig;

    % ??
    B = V * diag(d) * V';

    % ????
    B = (B + B') / 2;

end