%% experiment II: classical conditioning
% description goes here

%% initialization
Jeong2022_preface;

%% used fixed random seed
rng(0);

%% key assumptions
use_clicks = 0;
use_cs_offset = 0;

%% experiment parameters
pre_cs_delay = 1;
cs_dur = 2;
trace_dur = 1;
iti_delay = 3;
iti_mu = 60;
iti_max = 190;

%% simulation parameters
n_trials_per_run = 30;
n_runs = 9;
n_trials = n_trials_per_run * n_runs;

%% conditioned stimuli (CS)
cs_plus_proportion = .5;
cs = categorical(rand(n_trials,1)<=cs_plus_proportion,[0,1],{'CS-','CS+'});
cs_set = categories(cs);
n_cs = numel(cs_set);
cs_plus_flags = cs == 'CS+';
cs_clrs = [.9,.1,.15;.05,.45,.75];

%% time

% trial time
trial_dur = pre_cs_delay + cs_dur + trace_dur + iti_delay;
max_trial_dur = trial_dur + iti_max;
trial_time = (0 : dt : max_trial_dur - dt) - pre_cs_delay;
n_states_per_trial = numel(trial_time);
trial_state_edges = linspace(0,max_trial_dur,n_states_per_trial+1);

% simulation time
dur = max_trial_dur * n_trials;
time = 0 : dt : dur - dt;
n_states = numel(time);
state_edges = linspace(0,dur,n_states+1);

%% inter-trial-intervals
iti_pd = truncate(makedist('exponential','mu',iti_mu),0,iti_max);
iti = random(iti_pd,n_trials,1);
iti = dt * round(iti / dt);

%% inter-trial-onset-intervals
itoi = trial_dur + iti;

%% trial onset times
trial_onset_times = cumsum(itoi);
trial_onset_times = dt * round(trial_onset_times / dt);

%% simulation time
% dur = trial_onset_times(end) + max_trial_dur;
% dur = (n_trials_per_run * max_trial_dur) * ...
%     ceil(dur / (n_trials_per_run * max_trial_dur));
% time = 0 : dt : dur - dt;
% n_states = numel(time);
% state_edges = linspace(0,dur,n_states+1);

%% CS onset times
cs_plus_onset_times = trial_onset_times(cs_plus_flags) + pre_cs_delay;
cs_minus_onset_times = trial_onset_times(~cs_plus_flags) + pre_cs_delay;
cs_plus_onset_counts = histcounts(cs_plus_onset_times,state_edges);
cs_minus_onset_counts = histcounts(cs_minus_onset_times,state_edges);

%%
% cs_plus_onset_times = cs_plus_onset_times + linspace(0,1,5) .* cs_dur;
% cs_plus_onset_times = cs_plus_onset_times(:);

%% CS offset times
cs_plus_offset_times = cs_plus_onset_times + cs_dur;
cs_minus_offset_times = cs_minus_onset_times + cs_dur;
cs_plus_offset_counts = histcounts(cs_plus_offset_times,state_edges);
cs_minus_offset_counts = histcounts(cs_minus_offset_times,state_edges);

%% click times
click_trial_times = nan(n_trials,1);
click_trial_times(cs_plus_flags) = pre_cs_delay + cs_dur + trace_dur;
click_counts_2d = histcounts2(1:n_trials,click_trial_times',...
    'xbinedges',1:n_trials+1,...
    'ybinedges',trial_state_edges);
click_times = click_trial_times + trial_onset_times;
click_times = dt * round(click_times / dt);
click_counts = histcounts(click_times,state_edges);
[~,click_state_idcs] = ...
    min(abs(time - click_times(cs_plus_flags)),[],2);
n_clicks = numel(click_state_idcs);

%% reaction times
reaction_times = repmat(.5,n_trials,1);
% reaction_times = linspace(1,.1,n_trials)';
reaction_times = max(reaction_times,dt*2);
if ~use_clicks
    reaction_times = 0;
end

%% reward times
reward_trial_times = click_trial_times + reaction_times;
reward_trial_times = dt * floor(reward_trial_times / dt);
reward_counts_2d = histcounts2(1:n_trials,reward_trial_times',...
    'xbinedges',1:n_trials+1,...
    'ybinedges',trial_state_edges);
reward_times = click_times + reaction_times;
reward_times = dt * round(reward_times / dt);
reward_counts = histcounts(reward_times,state_edges);
[~,reward_state_idcs] = ...
    min(abs(time - reward_times(cs_plus_flags)),[],2);
n_rewards = numel(reward_state_idcs);

%% microstimuli & eligibility traces

% microstimuli
stimulus_trace = stimulustracefun(y0,tau,time)';
mus = linspace(1,0,n);
microstimuli = microstimulusfun(stimulus_trace,mus,sigma);

%% UNCOMMENT TO REPLACE MICROSTIMULI WITH COMPLETE SERIAL COMPOUND
% csc = zeros(n_states,n);
% pulse_duration = .5;
% pulse_length = floor(pulse_duration / dt);
% for ii = 1 : n
%     idcs = (1 : pulse_length) + (ii - 1) * pulse_length;
%     csc(idcs,ii) = 1;
% end
% microstimuli = csc;
% microstimuli = microstimuli / max(sum(microstimuli,2));

%% concatenate all stimulus times
stimulus_times = {...
    cs_plus_onset_times,...
    cs_minus_onset_times,...
    reward_times};
if use_clicks
    stimulus_times = [...
        stimulus_times,...
        click_times];
end
if use_cs_offset
    stimulus_times = [...
        stimulus_times,...
        cs_plus_offset_times,...
        cs_minus_offset_times];
end

%% TD learning
[state,value,rpe,exp2_weights] = tdlambda(...
    time,stimulus_times,reward_times,microstimuli,[],...
    'alpha',alpha,...
    'gamma',gamma,...
    'lambda',lambda);

%% compute 'DA signal'
padded_rpe = padarray(rpe,dlight_kernel.nbins/2,0);
da = conv(padded_rpe(1:end-1),dlight_kernel.pdf,'valid');

%% reshape from trial-less time series to STATES x TRIALS x RUNS tensors

% preallocation
stimulus_tensor = nan(n_states_per_trial,n_trials_per_run,n_runs);
value_tensor = nan(n_states_per_trial,n_trials_per_run,n_runs);
rpe_tensor = nan(n_states_per_trial,n_trials_per_run,n_runs);
da_tensor = nan(n_states_per_trial,n_trials_per_run,n_runs);

% iterate through runs
for ii = 1 : n_runs
    
    % iterate through trials
    for jj = 1 : n_trials_per_run
        trial_idx = jj + (ii - 1) * n_trials_per_run;
        onset_idx = find(time >= trial_onset_times(trial_idx),1);
        if trial_idx < n_trials
            offset_idx = find(time >= trial_onset_times(trial_idx+1),1);
        else
            offset_idx = find(time >= trial_onset_times(trial_idx) + trial_dur,1);
        end
        idcs = onset_idx : offset_idx - 1;
        n_idcs = numel(idcs);
        stimulus_tensor(1:n_idcs,jj,ii) = ...
            cs_plus_onset_counts(idcs) + ...
            cs_minus_onset_counts(idcs) + ...
            cs_plus_offset_counts(idcs) * use_cs_offset + ...
            cs_minus_offset_counts(idcs) * use_cs_offset + ...
            click_counts(idcs) * use_clicks + ...
            reward_counts(idcs);
        value_tensor(1:n_idcs,jj,ii) = value(idcs);
        rpe_tensor(1:n_idcs,jj,ii) = rpe(idcs);
        da_tensor(1:n_idcs,jj,ii) = da(idcs);
    end
end

%% figure 2: experiment II

% figure initialization
figure(figopt,...
    'name','experiment II: test 3');

% axes initialization
n_rows = 3 + 2;
n_cols = n_runs;
sp_stimulus = gobjects(1,n_runs);
sp_da_mu = gobjects(1,n_runs);
sp_da = gobjects(1,n_runs);
sp_value_mu = gobjects(1,n_runs);
sp_value = gobjects(1,n_runs);
for ii = 1 : n_runs
    sp_stimulus(ii) = subplot(n_rows,n_cols,ii+n_cols*0);
    sp_da_mu(ii) = subplot(n_rows,n_cols,ii+n_cols*1);
    sp_da(ii) = subplot(n_rows,n_cols,ii+n_cols*2);
    sp_value_mu(ii) = subplot(n_rows,n_cols,ii+n_cols*3);
    sp_value(ii) = subplot(n_rows,n_cols,ii+n_cols*4);
end

% concatenate axes
sps_stages = [...
    sp_stimulus;...
    sp_da_mu;...
    sp_da;...
    sp_value_mu;...
    sp_value;...
    ];
sps = [...
    sps_stages(:)',...
    ];

% axes settings
set(sps,axesopt);
set(sps_stages,...
    'xlim',[-pre_cs_delay,trial_dur+iti_delay]);
set([sp_da,sp_value],...
    'colormap',bone(2^8));

% axes titles
title(sp_stimulus(1),'Early in training');
title(sp_stimulus(end),'Late in training');

% axes labels
arrayfun(@(ax)xlabel(ax,'Time (s)'),sp_stimulus);
arrayfun(@(ax)ylabel(ax,'Stimulus trace (a.u.)'),sp_stimulus);
arrayfun(@(ax)xlabel(ax,'Time (s)'),[sp_da_mu;sp_da]);
arrayfun(@(ax)ylabel(ax,'DA (a.u.)'),[sp_da_mu;sp_da]);
arrayfun(@(ax)xlabel(ax,'Time (s)'),[sp_value_mu;sp_value]);
arrayfun(@(ax)ylabel(ax,'Value (a.u.)'),[sp_value_mu;sp_value]);

% iterate through runs
for ii = 1 : n_runs
    run_idcs = (1 : n_trials_per_run) + (ii - 1) * n_trials_per_run;
    [~,sorted_idcs] = sortrows(cs(run_idcs));
    
    % plot stimulus trace
    clim = quantile(stimulus_tensor,[0,1],'all')';
    imagesc(sp_stimulus(ii),...
        [trial_time(1),trial_time(end)],[],...
        stimulus_tensor(:,sorted_idcs,ii)',clim);
    
    % plot DA trace
    clim = quantile(da_tensor,[0,1],'all')';
    imagesc(sp_da(ii),...
        [trial_time(1),trial_time(end)],[],...
        da_tensor(:,sorted_idcs,ii)',clim);
    
    % plot value trace
    clim = quantile(value_tensor,[0,1],'all')';
    imagesc(sp_value(ii),...
        [trial_time(1),trial_time(end)],[],...
        value_tensor(:,sorted_idcs,ii)',clim);
    
    % iterate through CS conditions
    for jj = 1 : n_cs
        cs_flags = cs(run_idcs) == cs_set{jj};
        
        % compute & plot average DA conditioned on CS
        da_mu = nanmean(da_tensor(:,cs_flags,ii),2);
        plot(sp_da_mu(ii),trial_time,da_mu,...
            'color',cs_clrs(jj,:),...
            'linewidth',1);
%         rpe_mu = nanmean(rpe_tensor(:,cs_flags,ii),2);
%         stem(sp_da_mu(ii),trial_time,rpe_mu,...
%             'marker','none',...
%             'color',cs_clrs(jj,:));
        
        % compute & plot average value conditioned on CS
        value_mu = nanmean(value_tensor(:,cs_flags,ii),2);
        plot(sp_value_mu(ii),trial_time,value_mu,...
            'color',cs_clrs(jj,:),...
            'linewidth',1);
    end
end

% legend
[leg,icons] = legend(sp_da_mu(1),cs_set,...
    'location','northeast',...
    'box','off',...
    'autoupdate','off');
icons(3).XData = icons(3).XData + [1,0] * .3;
icons(5).XData = icons(5).XData + [1,0] * .3;

% axes linkage
arrayfun(@(ax1,ax2,ax3,ax4,ax5)linkaxes([ax1,ax2,ax3,ax4,ax5],'x'),...
    sp_stimulus,sp_da_mu,sp_da,sp_value_mu,sp_value);
linkaxes(sp_da_mu,'y');
linkaxes(sp_value_mu,'y');

% iterate through runs
for ii = 1 : n_runs

    % plot CS onset
    plot(sp_da_mu(ii),[0,0],ylim(sp_da_mu(ii)),...
        'color','k',...
        'linestyle','--',...
        'linewidth',1);
    plot(sp_value_mu(ii),[0,0],ylim(sp_value_mu(ii)),...
        'color','k',...
        'linestyle','--',...
        'linewidth',1);
    
    % plot CS offset
    if use_cs_offset
        plot(sp_da_mu(ii),[0,0]+cs_dur,ylim(sp_da_mu(ii)),...
            'color','k',...
            'linestyle','--',...
            'linewidth',1);
        plot(sp_value_mu(ii),[0,0]+cs_dur,ylim(sp_value_mu(ii)),...
            'color','k',...
            'linestyle','--',...
            'linewidth',1);
    end
    
    % plot offset of the trace period
    plot(sp_da_mu(ii),[0,0]+cs_dur+trace_dur,ylim(sp_da_mu(ii)),...
        'color','k',...
        'linestyle','--',...
        'linewidth',1);
    plot(sp_value_mu(ii),[0,0]+cs_dur+trace_dur,ylim(sp_value_mu(ii)),...
        'color','k',...
        'linestyle','--',...
        'linewidth',1);
end

% annotate model parameters
annotateModelParameters;