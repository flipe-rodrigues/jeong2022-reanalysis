%% save settings
want2savedata = 1;
want2savepanels = 1;

%% normalization settings
want2renormalize = 1;

%% directory settings
mfilepath = mfilename('fullpath');
if contains(mfilepath,'LiveEditorEvaluationHelper')
    mfilepath = matlab.desktop.editor.getActiveFilename;
end
root_path = fileparts(fileparts(fileparts(mfilepath)));
panel_path = fullfile(root_path,'panels');
data_path = fullfile(root_path,'data');
experiments_path = fullfile(data_path,'experiments');
mice_path = fullfile(data_path,'mice');
mice_dir = dir(mice_path);
mice_dir = mice_dir(cellfun(@(x)~contains(x,'.'),{mice_dir.name}));

%% mouse settings
mouse_ids = {mice_dir.name}';
mouse_ids = mouse_ids([2,4,6,7,8,1,3,5]);
n_mice = numel(mouse_ids);

%% acquisition settings
fs = 120;
dt = 1 / fs;

%% smoothing kernels
lickrate_kernel = gammakernel(...
    'peakx',.15,...
    'binwidth',dt);