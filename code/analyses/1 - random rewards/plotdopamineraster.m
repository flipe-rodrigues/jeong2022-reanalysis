function plotdopamineraster(data,meta,opt)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
da_mat = data.da.(opt.alignment);
imagesc(meta.epochs.(opt.alignment),[0,size(da_mat,1)]+[1,-1]*.5,...
    da_mat,quantile(da_mat,[.001,.999],'all')');
end

