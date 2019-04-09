% SUMMARY:  Calculate ln p(xn|zn) of Gaussian
% AUTHOR:   QIUQIANG KONG
% Created:  17-11-2015
% Modified: 25-11-2015 Add annotation
% -----------------------------------------------------------
% input:
%   X       size: N*p
%   phi
%     mu      size: p*M
%     Sigma   size: p*p*M
% output:
%   ln p(xn|zn)
% ===========================================================
function logp_xn_given_zn = Gauss_logp_xn_given_zn(X, mu,sigma)
[N,p] = size(X);
[p,M] = size(mu);
logp_xn_given_zn = zeros(N,M);
for m = 1:M
    x_minus_mu = bsxfun(@minus, X, mu(:,m)');
    logp_xn_given_zn(:,m) = -0.5*p*log(2*pi) - 0.5*log(det(sigma(:,:,m))) - 0.5 * sum(x_minus_mu * inv(sigma(:,:,m)) .* x_minus_mu, 2);
end
end
