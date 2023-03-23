% converts from/to temporal discount factor to/from time constant [s]
gammafun = @(dt,eta) exp(-dt/eta);
etafun = @(dt,gamma) -dt/log(gamma);

%% microstimuli definition
stimulustracefun = @(y0,tau,t) ...
    y0 * tau .^ t;
gaussianbasisfun = @(y,mu,sigma) ...
    1 / sqrt(2 * pi) * exp(-(y - mu).^2 / (2 * sigma^2));
microstimulusfun = @(yt,mu,sigma) ...
    gaussianbasisfun(yt,mu,sigma) .* yt;