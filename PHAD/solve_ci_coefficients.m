function [cG, cU, info] = solve_ci_coefficients(H, dg, m, C)
% solve_ci_coefficients
%
% Computes cG and cA without explicit matrix inverse.
%
% Case 1: MQMDPA
%   [cG, cA] = solve_ci_coefficients(H, g, B)
%
%   H \ cG = -dg
%   H \ cA = m
%
%   Usually:
%       m = Aeq.'
%   where Aeq is Nc x Npa.
%
% Case 2: CEF
%   [cG, cA] = solve_ci_coefficients(H, g, B, C)
%
%   [H C.'; C 0] [c; eta] = RHS
%
% Inputs:
%   H : Npa x Npa Hessian
%   dg : Npa x 1 gradient
%   m : Npa x Nc RHS matrix, e.g. dMdy or Aeq.'
%   C : optional constraint matrix, ns x ny
%
% Outputs:
%   cG   : Npa x 1
%   cA   : NPa x Nc
%   info : diagnostic info

    dg = dg(:);
    ny = length(dg);
    

    if size(H,1) ~= ny || size(H,2) ~= ny
        error('H must be ny x ny, where ny = length(g).');
    end

    if size(m,1) ~= ny
        error('m must have ny rows.');
    end

    % Multiple RHS: first column for cG, remaining columns for cA
    RHS0 = [-dg, m];

    if nargin < 4 || isempty(C)
        % MQMDPA case: no internal normalization constraint
        SOL = H \ RHS0;

        cG = SOL(:,1);
        cU = SOL(:,2:end);

        info.mode = 'unconstrained';
        info.matrix_size = size(H);
        info.cond_est = cond(H);

    else
        % CEF case: constrained site-fraction variables
        ns = size(C,1);

        if size(C,2) ~= ny
            error('C must have ny columns.');
        end

        P = [H, C.';
             C, zeros(ns,ns)];

        RHS = [RHS0;
               zeros(ns, size(RHS0,2))];

        SOL = P \ RHS;

        cG = SOL(1:ny,1);
        cU = SOL(1:ny,2:end);

        info.mode = 'constrained';
        info.matrix_size = size(P);
        info.cond_est = cond(P);
    end
end