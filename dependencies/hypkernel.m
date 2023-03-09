function k = hypkernel(varargin)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
     
    p = inputParser;
    p.addParameter('med',10);
    p.addParameter('binwidth',2);
    p.parse(varargin{:});   
    
    k = p.Results;
    k.type = 'hyperbolic';
    k.k = 21.365 / k.med;
    k.nbins = (1 / k.k) / k.binwidth * 1e3;
    k.paddx = [-1,1] / 2 * k.nbins * k.binwidth;
    k.bins = -k.nbins / 2 + 1 : k.nbins / 2;
    k.x = k.bins * k.binwidth;
    k.pdf = 1 ./ (1 + k.k .* k.x);
    k.pdf(k.x < 0) = 0;
    k.pdf = k.pdf / sum(k.pdf);
end
