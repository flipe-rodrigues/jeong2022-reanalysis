%% initialization
close all;
clear;
clc;

%% directory settings
root_path = fileparts(pwd);
data_path = fullfile(root_path,'data');
data_dir = dir(data_path);
data_dir = data_dir(cellfun(@(x)~contains(x,'.'),{data_dir.name}));
save_path = fullfile(root_path,'figures');
want2save = false;

%% normalization settings
want2renormalize = true;

%% mouse settings
mouse_ids = {data_dir.name};
mouse_ids = mouse_ids([6,2,3,4,5,1,7,8]);
n_mice = numel(mouse_ids);

%% acquisition settings
fs = 120;
dt = 1 / fs;

%% smoothing kernels
lickrate_kernel = gammakernel(...
    'peakx',.15,...
    'binwidth',dt);
