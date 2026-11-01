function [G, S, H, Cp] = G0_pure(Par, comp, phase, T)

    if ~isfield(Par, 'G0')
        error('G0_pure:MissingG0','Par.G0 is missing.');
    end

    if ~isfield(Par.G0, comp)
        error('G0_pure:MissingComp','Missing G0 for component/species: %s', comp);
    end

    if ~isfield(Par.G0.(comp), phase)
        error('G0_pure:MissingPhase','Missing phase %s for %s', phase, comp);
    end

    segs = Par.G0.(comp).(phase);

    G = nan(size(T));
    S = nan(size(T));
    H = nan(size(T));
    dSdT = nan(size(T));
    Cp = nan(size(T));

    for k = 1:numel(segs)
        s = segs(k);

        idT = (T >= s.Tmin) & (T <= s.Tmax);

        if ~any(idT)
            continue;
        end

        Tk = T(idT);

        G(idT) = s.a ...
               + s.b .* Tk ...
               + s.c .* Tk .* log(Tk) ...
               + s.d .* Tk.^2 ...
               + s.e .* Tk.^3 ...
               + s.f ./ Tk ...
               + s.g .* Tk.^1.5 ...
               + s.h .* Tk.^0.5 ...
               + s.i .* log(Tk);
           
        S(idT) = -(s.b + s.c ...
               + s.c .* log(Tk) ...
               + 2 * s.d .* Tk ...
               + 3 * s.e .* Tk.^2 ...
               - s.f ./ Tk.^2 ...
               + 1.5 * s.g .* Tk.^0.5 ...
               + 0.5 * s.h .* Tk.^(-0.5) ...
               + s.i .* Tk.^(-1));
           
        dSdT(idT) = -(s.c ./ Tk ... 
             + 2 .* s.d ...
             + 6 .* s.e .* Tk ...
             + 2 .* s.f ./ Tk.^3 ...
             + 0.75 .* s.g ./ Tk.^0.5 ...
             - 0.25 .* s.h ./ Tk.^1.5 ...
             - s.i ./ Tk.^2); 
           
        H(idT) = G(idT) + Tk .* S(idT);
        
        Cp(idT) = Tk .* dSdT(idT);
        
    end

    if any(isnan(G))
        error('G0_pure:TemperatureOutOfRange', ...
            'Temperature outside valid range for %s %s.', comp, phase);
    end
    
    if any(Cp < 0)
        warning('G0_pure:NegativeCp','Negative Cp detected for %s %s.', comp, phase);
    end
end