% I need a model that has different rhos for gain trials than loss trials.
% I don't think that this implementation of prospecet theory gives me that.

c1 = 10;
c2 = 5;
p1 = 0.5;
rho = 0;
alpha = 1;%.7; % alpha is the weight of the probability

u_CARA   = @(c, rho)   (rho ~= 0) .* (1 - exp(-rho .* c)) ./ (   rho + (rho == 0) .* 0.000001)  + ...
                        (rho == 0) .*                  c;

u = u_CARA;

EU_CARA   = @(c1, c2, p1, rho)  p1 .* u_CARA(c1, rho) + (1 - p1) .* u_CARA(c2, rho);

EU_Pr_CARA  =  @(c1, c2, p1, rho, alpha)  (1 - w_Pr(1 - p1, alpha))    .* u_CARA(c1, rho) + ...
                                               w_Pr(1 - p1, alpha)     .* u_CARA(c2, rho);
% Probability weighting as proposed by Tversky & Kahneman (J. Risk & Uncertainty, 1992)
 w_TK    =  @(p, alpha)                (p > 0) .* (p.^alpha ./ (p.^alpha + (1 - p).^alpha).^(1./alpha));
% EU according to cumulative prospect theory:
EU_TK    =  @(c1, c2, p1, rho, alpha)  (1 - w_TK(1 - p1, alpha))       .* u(c1, rho) + ...
                                            w_TK(1 - p1, alpha)        .* u(c2, rho);

eu = EU_TK(c1, c2, p1, rho, alpha)