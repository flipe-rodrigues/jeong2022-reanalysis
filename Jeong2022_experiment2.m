%% experiment II: classical conditioning
% description goes here

%% initialization
Jeong2022_preface;

%% used fixed random seed
rng(0);

%% key assumptions
use_clicks = 0;
use_cs_offset = 1;

%% experiment parameters
pre_cs_delay = 1;
cs_dur_set = 2;
n_cs_durs = numel(cs_dur_set);
trace_dur = 1;
iti_delay = 3;
iti_mu = 30;
iti_max = 90;

%% analysis parameters
baseline_period = [-1,0];
early_period = [0,1];
late_period = [-1,0];

%% simulation parameters
n_trials_per_run = 100;
n_runs = 9;
n_trials = n_trials_per_run * n_runs;

%% conditioned stimuli (CS)
cs_dur = repmat(cs_dur_set(1),n_trials,1);
cs_plus_proportion = .5;
cs = categorical(rand(n_trials,1)<=cs_plus_proportion,[0,1],{'CS-','CS+'});
cs_set = categories(cs);
n_cs = numel(cs_set);
cs_plus_flags = cs == 'CS+';

%% trial time
trial_dur = pre_cs_delay + cs_dur + trace_dur + iti_delay;
max_trial_dur = max(trial_dur) + iti_max;
trial_time = (0 : dt : max_trial_dur - dt) - pre_cs_delay;
n_states_per_trial = numel(trial_time);
trial_state_edges = linspace(0,max_trial_dur,n_states_per_trial+1);

%% inter-trial-intervals
iti_pd = truncate(makedist('exponential','mu',iti_mu),0,iti_max);
iti = random(iti_pd,n_trials,1);
iti = dt * round(iti / dt);
n_bins = round(max(iti) / iti_mu) * 10;
iti_edges = linspace(0,max(iti),n_bins);
iti_counts = histcounts(iti,iti_edges);
iti_counts = iti_counts ./ nansum(iti_counts);
iti_pdf = pdf(iti_pd,iti_edges);
iti_pdf = iti_pdf ./ nansum(iti_pdf);

%% inter-trial-onset-intervals
itoi = trial_dur + iti;

%% trial onset times
trial_onset_times = cumsum(itoi);
trial_onset_times = dt * round(trial_onset_times / dt);

%% simulation time
dur = trial_onset_times(end) + max_trial_dur;
time = 0 : dt : dur - dt;
n_states = numel(time);
state_edges = linspace(0,dur,n_states+1);

%% CS onset times
cs_plus_onset_times = trial_onset_times(cs_plus_flags) + pre_cs_delay;
cs_minus_onset_times = trial_onset_times(~cs_plus_flags) + pre_cs_delay;
cs_plus_onset_counts = histcounts(cs_plus_onset_times,state_edges);
cs_minus_onset_counts = histcounts(cs_minus_onset_times,state_edges);

%% CS offset times
cs_plus_offset_times = cs_plus_onset_times + cs_dur(cs_plus_flags);
cs_minus_offset_times = cs_minus_onset_times + cs_dur(~cs_plus_flags);
cs_plus_offset_counts = histcounts(cs_plus_offset_times,state_edges);
cs_minus_offset_counts = histcounts(cs_minus_offset_times,state_edges);

%% CS flags
cs_plus_ison = sum(...
    time >= cs_plus_onset_times & ...
    time <= cs_plus_offset_times,1)';
cs_minus_ison = sum(...
    time >= cs_minus_onset_times & ...
    time <= cs_minus_offset_times,1)';

%% click times
click_trial_times = nan(n_trials,1);
click_trial_times(cs_plus_flags) = pre_cs_delay + cs_dur(cs_plus_flags) + trace_dur;
click_times = click_trial_times + trial_onset_times;
click_times = dt * round(click_times / dt);
click_counts = histcounts(click_times,state_edges);
n_clicks = sum(click_counts);

%% reaction times
reaction_times = repmat(.5,n_trials,1);
% reaction_times = linspace(1,.1,n_trials)';
reaction_times = max(reaction_times,dt*2);
if ~use_clicks
    reaction_times = 0;
end

%% reward times
reward_times = click_times + reaction_times;
reward_times = dt * round(reward_times / dt);
reward_counts = histcounts(reward_times,state_edges);
n_rewards = sum(reward_counts);

%% microstimuli
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
    reward_times,...
    cs_plus_onset_times,...
    cs_minus_onset_times};
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

%% concatenate all state flags
cs_flags = [...
    cs_plus_ison,...
    cs_minus_ison];

%% TD learning
[state,value,rpe,exp2_weights] = tdlambda(...
    time,[],stimulus_times,reward_times,microstimuli,[],...
    'alpha',alpha,...
    'gamma',gamma,...
    'lambda',lambda);

%% compute 'DA signal'
padded_rpe = padarray(rpe,dlight_kernel.nbins/2,0);
da = conv(padded_rpe(1:end-1),dlight_kernel.pdf,'valid');

%% get reward-aligned snippets of DA & value signals
[da_baseline_snippets,da_baseline_time] = ...
    signal2eventsnippets(time,da,cs_plus_onset_times,baseline_period,dt);
[da_early_snippets,da_early_time] = ...
    signal2eventsnippets(time,da,cs_plus_onset_times,early_period,dt);
[da_late_snippets,da_late_time] = ...
    signal2eventsnippets(time,da,reward_times(cs_plus_flags),late_period,dt);

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
            offset_idx = find(time >= trial_onset_times(trial_idx) + trial_dur(trial_idx),1);
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
n_rows = 5;
n_cols = n_runs;
sp_cstype = subplot(n_rows,n_cols,1);
sp_state = subplot(n_rows,n_cols,2:n_cols-2);
sp_iti = subplot(n_rows,n_cols,n_cols-1:n_cols);
sp_da_mu = gobjects(1,n_runs);
sp_da = gobjects(1,n_runs);
sp_value_mu = gobjects(1,n_runs);
sp_value = gobjects(1,n_runs);
for ii = 1 : n_runs
    sp_da_mu(ii) = subplot(n_rows,n_cols,ii+n_cols*1);
    sp_da(ii) = subplot(n_rows,n_cols,ii+n_cols*2);
    sp_value_mu(ii) = subplot(n_rows,n_cols,ii+n_cols*3);
    sp_value(ii) = subplot(n_rows,n_cols,ii+n_cols*4);
end

% concatenate axes
sps_stages = [...
    sp_da_mu;...
    sp_da;...
    sp_value_mu;...
    sp_value;...
    ];
sps = [...
    sp_cstype;...
    sp_state;...
    sp_iti;...
    sps_stages(:);...
    ];

% axes settings
arrayfun(@(ax)set(ax.XAxis,'exponent',0),sps);
set(sps,axesopt);
set(sps_stages,...
    'xlim',[-pre_cs_delay,unique(trial_dur)+iti_delay]);
set(sp_cstype,...
    'xlim',[0,1]+[-1,1],...
    'xtick',[0,1],...
    'xticklabel',cs_set);
set(sp_state,...
    'ticklength',axesopt.ticklength/n_cols);
set([sp_state,sp_iti],...
    'ytick',[]);

% axes titles
title(sp_da_mu(1),'Early in training');
title(sp_da_mu(end),'Late in training');

% axes labels
ylabel(sp_cstype,'Count');
xlabel(sp_state,'Time (s)');
ylabel(sp_state,'State feature #');
xlabel(sp_iti,'ITI (s)');
ylabel(sp_iti,'PDF');
arrayfun(@(ax)xlabel(ax,'Time (s)'),[sp_da_mu;sp_da]);
arrayfun(@(ax)ylabel(ax,'DA (a.u.)'),[sp_da_mu;sp_da]);
arrayfun(@(ax)xlabel(ax,'Time (s)'),[sp_value_mu;sp_value]);
arrayfun(@(ax)ylabel(ax,'Value (a.u.)'),[sp_value_mu;sp_value]);

% plot CS distribution
patch(sp_cstype,...
    0+[-1,1,1,-1]*1/4,[0,0,1,1]*sum(~cs_plus_flags),cs_minus_clr,...
    'edgecolor',cs_minus_clr,...
    'linewidth',1.5,...
    'facealpha',2/3);
patch(sp_cstype,...
    1+[-1,1,1,-1]*1/4,[0,0,1,1]*sum(cs_plus_flags),cs_plus_clr,...
    'edgecolor',cs_plus_clr,...
    'linewidth',1.5,...
    'facealpha',2/3);

% plot state features
n_trials2plot = 15;
time_flags = ...
    time >= trial_onset_times(1) & ...
    time < trial_onset_times(n_trials2plot + 1);
imagesc(sp_state,...
    time(time_flags)+dt/2,[],state(time_flags,:)');

% plot ITI distribution
stem(sp_iti,iti_mu,max([iti_counts,iti_pdf]),...
    'color','k',...
    'marker','v',...
    'markersize',10,...
    'markerfacecolor','k',...
    'markeredgecolor','none',...
    'linewidth',2);
histogram(sp_iti,...
    'binedges',iti_edges,...
    'bincounts',iti_counts,...
    'facecolor','w',...
    'edgecolor','k',...
    'facealpha',1,...
    'linewidth',1);
plot(sp_iti,iti_edges,iti_pdf,...
    'color','k',...
    'linewidth',2);

% iterate through runs
for ii = 1 : n_runs
    run_idcs = (1 : n_trials_per_run) + (ii - 1) * n_trials_per_run;
    [~,sorted_idcs] = sortrows(cs(run_idcs));
    
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
    
    % compute & plot average DA conditioned on CS
    da_minus_mu = nanmean(da_tensor(:,~cs_plus_flags(run_idcs),ii),2);
    plot(sp_da_mu(ii),trial_time,da_minus_mu,...
        'color',cs_minus_clr,...
        'linewidth',1);
    da_plus_mu = nanmean(da_tensor(:,cs_plus_flags(run_idcs),ii),2);
    plot(sp_da_mu(ii),trial_time,da_plus_mu,...
        'color',cs_plus_clr,...
        'linewidth',1);
    %         rpe_mu = nanmean(rpe_tensor(:,cs_flags,ii),2);
    %         plot(sp_da_mu(ii),trial_time,rpe_mu,...
    %             'color',cs_clrs(jj,:),...
    %             'linewidth',1);
    
    % compute & plot average value conditioned on CS
    value_minus_mu = nanmean(value_tensor(:,~cs_plus_flags(run_idcs),ii),2);
    plot(sp_value_mu(ii),trial_time,value_minus_mu,...
        'color',cs_minus_clr,...
        'linewidth',1);
    value_plus_mu = nanmean(value_tensor(:,cs_plus_flags(run_idcs),ii),2);
    plot(sp_value_mu(ii),trial_time,value_plus_mu,...
        'color',cs_plus_clr,...
        'linewidth',1);
end

% legend
[leg,icons] = legend(sp_da_mu(1),cs_set,...
    'location','northeast',...
    'box','off',...
    'autoupdate','off');
icons(3).XData = icons(3).XData + [1,0] * .3;
icons(5).XData = icons(5).XData + [1,0] * .3;

% axes linkage
arrayfun(@(ax1,ax2,ax3,ax4)linkaxes([ax1,ax2,ax3,ax4],'x'),...
    sp_da_mu,sp_da,sp_value_mu,sp_value);
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
        plot(sp_da_mu(ii),[0,0]+unique(cs_dur),ylim(sp_da_mu(ii)),...
            'color','k',...
            'linestyle','--',...
            'linewidth',1);
        plot(sp_value_mu(ii),[0,0]+unique(cs_dur),ylim(sp_value_mu(ii)),...
            'color','k',...
            'linestyle','--',...
            'linewidth',1);
    end
    
    % plot offset of the trace period
    plot(sp_da_mu(ii),[0,0]+unique(cs_dur)+trace_dur,ylim(sp_da_mu(ii)),...
        'color','k',...
        'linestyle','--',...
        'linewidth',1);
    plot(sp_value_mu(ii),[0,0]+unique(cs_dur)+trace_dur,ylim(sp_value_mu(ii)),...
        'color','k',...
        'linestyle','--',...
        'linewidth',1);
end

% annotate model parameters
annotateModelParameters;