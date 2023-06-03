function k = expkernel(varargin)
    %EXPKERNEL Creates an exponential kernel
    %   K = EXPKERNEL('mus',MUS,'binwidth',BINWIDTH) creates an exponential
    %   kernel with a mean of MUS ms and a resolution of BINWIDTH ms.
     
    p = inputParser;
    p.addParameter('mus',10);
    p.addParameter('binwidth',2e-3);
    p.parse(varargin{:});   
    
    k = p.Results;
    k.type = 'exponential';
    k.n = numel(k.mus);
    k.taus = 1 ./ k.mus;
    k.nbins = max(k.mus) / k.binwidth * 20;
    k.paddx = [-1,1] / 2 * k.nbins * k.binwidth;
    k.bins = -k.nbins / 2 + 1 : k.nbins / 2;
    k.x = k.bins * k.binwidth;
    for ii = 1 : k.n
        k.pdfs(ii,:) = exppdf(k.x,k.mus(ii));
        k.pdfs(ii,:) = k.pdfs(ii,:) / sum(k.pdfs(ii,:));
    end
    k.pdf = sum(k.pdfs,1);
    k.pdf = k.pdf / sum(k.pdf);
end
