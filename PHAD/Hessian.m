function H = Hessian(gradfun, x, h0)

    x = real(x(:));
    n = length(x);
    H = zeros(n,n);

    for k = 1:n

        hk = h0 * max(1, abs(x(k)));

        % ?? x(k)-hk ??
        if x(k) - hk <= 1e-16
            xp = x;
            x0 = x;

            xp(k) = x(k) + hk;

            gp = gradfun(xp);
            g0 = gradfun(x0);

            H(:,k) = (gp(:) - g0(:)) / hk;
        else
            xp = x;
            xm = x;

            xp(k) = x(k) + hk;
            xm(k) = x(k) - hk;

            gp = gradfun(xp);
            gm = gradfun(xm);

            H(:,k) = (gp(:) - gm(:)) / (2*hk);
        end
    end

    H = 0.5 * (H + H.');
end

% function H = Hessian(gradfun, x, h)
%     x = x(:);
%     n = length(x);
%     H = zeros(n,n);
% 
%     for j = 1:n
%         hj = h;
% 
%         xp = x;
%         xm = x;
% 
%         xp(j) = xp(j) + hj;
%         xm(j) = max(xm(j) - hj, eps);
% 
%         gp = gradfun(xp);
%         gm = gradfun(xm);
% 
%         H(:,j) = (gp(:) - gm(:)) / (2*hj);
%     end
% 
%     H = 0.5 * (H + H.');
% end

% function H = numericalHessian_from_energy(Efun, x, h)
% 
%     if nargin < 3
%         h = 1e-5;
%     end
% 
%     x = x(:);
%     n = length(x);
%     H = zeros(n,n);
% 
%     for i = 1:n
%         hi = h * max(abs(x(i)), 1);
%         for j = i:n
%             hj = h * max(abs(x(j)), 1);
% 
%             xp = x; xp(i) = xp(i) + hi; xp(j) = xp(j) + hj;
%             xm = x; xm(i) = xm(i) - hi; xm(j) = xm(j) - hj;
%             xpm = x; xpm(i) = xpm(i) + hi; xpm(j) = xpm(j) - hj;
%             xmp = x; xmp(i) = xmp(i) - hi; xmp(j) = xmp(j) + hj;
% 
%             fpp = Efun(xp);
%             fmm = Efun(xm);
%             fpm = Efun(xpm);
%             fmp = Efun(xmp);
% 
%             H(i,j) = (fpp + fmm - fpm - fmp) / (4*hi*hj);
%             H(j,i) = H(i,j);
%         end
%     end
% end