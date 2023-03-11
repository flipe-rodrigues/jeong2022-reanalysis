%% experiment I: unpredicted rewards
% description goes here

%% initialization
Jeong2022_preface;

%% used fixed random seed
rng(0);

%% key assumptions
use_clicks = 1;

%% analysis parameters
iri_cutoff = 3;
reward_period = [-.5,1];
baseline_period = [-2,-.5];

%% experiment parameters
iri_mu = 12;

%% simulation parameters
n_rewards = 500;

%% training stage settings
n_stages = 3;
stage_dur = 60 * 2;
early_clr = [0,0,0];
late_clr = [1,1,1] * .85;
stage_clrs = colorlerp([early_clr; late_clr],n_stages);

%% inter-reward-intervals
iri_pd = makedist('exponential','mu',iri_mu);
iri = random(iri_pd,n_rewards,1);

%% time
dur = sum(iri) + max(iri);
time = 0 : dt : dur - dt;
n_states = numel(time);
state_edges = linspace(0,dur,n_states+1);

%% click times
click_times = cumsum(iri);
click_times = unique(dt * round(click_times / dt));
click_counts = histcounts(click_times,state_edges);

%% reaction times
reaction_times = repmat(.5,n_rewards,1);
reaction_times = linspace(1,.1,n_rewards)';
% reaction_times = normalize01(click_times .^ -dt) + .1;
% reaction_times = normalize01(1 ./ (1 + .005 .* click_times)) + .1;
% reaction_times = exprnd(.5,n_clicks,1);
% reaction_times = exprnd(linspace(1.5,1/3,n_clicks)',n_clicks,1);
reaction_times = max(reaction_times,dt*2);

%% reward times
reward_times = click_times + reaction_times;
reward_times = dt * round(reward_times / dt);
reward_counts = histcounts(reward_times,state_edges);

%% inter-reward-intervals
n_bins = round(max(iri) / iri_mu) * 10;
iri_edges = linspace(0,max(iri),n_bins);
iri_counts = histcounts(iri,iri_edges);
iri_counts = iri_counts ./ nansum(iri_counts);
iri_pdf = pdf(iri_pd,iri_edges);
iri_pdf = iri_pdf ./ nansum(iri_pdf);

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

%% elibility traces
eligibility = zeros(n_states,n);
for ss = 2 : n_states
    eligibility(ss,:) = ...
        gamma * lambda * eligibility(ss-1,:) + microstimuli(ss-1,:);
end

%% concatenate all stimulus times
stimulus_times = {...
    reward_times};
if use_clicks
    stimulus_times = [...
        stimulus_times,...
        click_times];
end

%% TD learning
[state,value,rpe,exp1_weights] = tdlambda(...
    time,[],stimulus_times,reward_times,microstimuli,[],...
    'alpha',alpha,...
    'gamma',gamma,...
    'lambda',lambda,...
    'tau',tau,...
    'theta',theta);

%% compute 'DA signal'
padded_rpe = padarray(rpe,dlight_kernel.nbins/2,0);
da = conv(padded_rpe(1:end-1),dlight_kernel.pdf,'valid');
da = da / max(dlight_kernel.pdf);

%% get reward-aligned snippets of DA signal
[da_reward_snippets,da_reward_time] = ...
    signal2eventsnippets(time,da,reward_times,reward_period,dt);
[da_baseline_snippets,da_baseline_time] = ...
    signal2eventsnippets(time,da,reward_times,baseline_period,dt);

%% IRI-based reward selection
reward_flags = iri >= iri_cutoff;
da_reward_snippets(~reward_flags,:) = nan;
da_baseline_snippets(~reward_flags,:) = nan;

%% compute 'DA response' metric

% preallocation
da_response = nan(n_rewards,1);

% iterate through rewards
for ii = 1 : n_rewards
    da_response(ii) = ...
        sum(da_reward_snippets(ii,:)) - sum(da_baseline_snippets(ii,:));
end

%% compute training stage indices
stage_dur = min(stage_dur,dur/n_stages);
n_states_per_stage = floor(n_states * min(stage_dur/dur,1/n_stages));
n_stage_divisors = floor(n_states / n_states_per_stage);
divisor_state_idcs = (0 : n_stage_divisors - 1) * n_states_per_stage;
stage_divisor_idcs = floor(linspace(1,n_stage_divisors,n_stages));
stage_state_idcs = (1 : n_states_per_stage)' + ...
    divisor_state_idcs(stage_divisor_idcs);

%% figure 1: experiment I

% figure initialization
figure(figopt,...
    'name','experiment I: tests 1 & 2');

% axes initialization
n_rows = 4;
n_cols = n_stages + 2;
sp_stimulus = gobjects(1,n_stages);
sp_state = gobjects(1,n_stages);
sp_rpe = gobjects(1,n_stages);
sp_value = gobjects(1,n_stages);
for ii = 1 : n_stages
    sp_stimulus(ii) = subplot(n_rows,n_cols,ii+n_cols*0);
    sp_state(ii) = subplot(n_rows,n_cols,ii+n_cols*1);
    sp_rpe(ii) = subplot(n_rows,n_cols,ii+n_cols*2);
    sp_value(ii) = subplot(n_rows,n_cols,ii+n_cols*3);
end
sp_iri = subplot(n_rows,n_cols,n_cols-1+n_cols*0);
sp_reaction = subplot(n_rows,n_cols,n_cols+n_cols*0);
sp_microstimulus = subplot(n_rows,n_cols,n_cols-1+n_cols*1);
sp_eligibility = subplot(n_rows,n_cols,n_cols+n_cols*1);
sp_baseline = subplot(n_rows,n_cols,n_cols-1+n_cols*2);
sp_reward = subplot(n_rows,n_cols,n_cols+n_cols*2);
sp_test1 = subplot(n_rows,n_cols,n_cols-1+n_cols*3);
sp_test2 = subplot(n_rows,n_cols,n_cols+n_cols*3);

% concatenate axes
sp_stage = [...
    sp_stimulus;...
    sp_state;...
    sp_rpe;...
    sp_value...
    ];
sps = [...
    sp_stage(:);...
    sp_iri;...
    sp_reaction;...
    sp_microstimulus;...
    sp_eligibility;...
    sp_baseline;...
    sp_reward;...
    sp_test1;...
    sp_test2;...
    ];

% axes settings
arrayfun(@(ax)set(ax.XAxis,'exponent',0),sps);
set(sps,axesopt);
set([sp_iri,sp_microstimulus,sp_eligibility],...
    'xlim',[0,40]);
set(sp_reward,...
    'ycolor','none');

% axes titles
title(sp_stimulus(1),'Early in training');
title(sp_stimulus(end),'Late in training');
title(sp_baseline,'Baseline period');
title(sp_reward,'Reward period');
title(sp_test1,'Test I');
title(sp_test2,'Test II');

% axes labels
arrayfun(@(ax)xlabel(ax,'Time (s)'),sp_stimulus);
arrayfun(@(ax)ylabel(ax,'Stimulus trace (a.u.)'),sp_stimulus);
arrayfun(@(ax)xlabel(ax,'Time (s)'),sp_state);
arrayfun(@(ax)ylabel(ax,'State feature #'),sp_state);
arrayfun(@(ax)xlabel(ax,'Time (s)'),sp_rpe);
arrayfun(@(ax)ylabel(ax,'RPE (a.u.)'),sp_rpe);
arrayfun(@(ax)xlabel(ax,'Time (s)'),sp_value);
arrayfun(@(ax)ylabel(ax,'Value (a.u.)'),sp_value);
xlabel(sp_iri,'IRI (s)');
ylabel(sp_iri,'PDF');
xlabel(sp_reaction,'Time (s)');
ylabel(sp_reaction,'Reaction time (s)');
xlabel(sp_microstimulus,'Time (s)');
ylabel(sp_microstimulus,'Microstimulus (a.u.)');
xlabel(sp_eligibility,'Time (s)');
ylabel(sp_eligibility,'Eligibility traces (a.u.)');
xlabel(sp_baseline,'Time since reward (s)');
ylabel(sp_baseline,'DA (a.u.)');
xlabel(sp_reward,'Time since reward (s)');
ylabel(sp_reward,'DA (a.u.)');
xlabel(sp_test1,'Reward #');
ylabel(sp_test1,'DA response (a.u.)');
xlabel(sp_test2,'Previous IRI (s)');
ylabel(sp_test2,'DA response (a.u.)');

% iterate through training stages
for ii = 1 : n_stages
    idcs = stage_state_idcs(:,ii);
    time_win = [time(idcs(1)),time(idcs(end))];
    
    % plot stimulus trace
    if use_clicks
        stem(sp_stimulus(ii),...
            time(idcs),click_counts(idcs),...
            'color','k',...
            'marker','none');
    end
    stem(sp_stimulus(ii),...
        time(idcs),reward_counts(idcs),...
        'color',reward_clr,...
        'marker','none');
    
    % plot state features
    imagesc(sp_state(ii),time(idcs)+dt/2,[],state(idcs,:)');
    
    % plot RPE
    stem(sp_rpe(ii),time(idcs),rpe(idcs),...
        'color','k',...
        'marker','none');
    
    % plot DA signal
    plot(sp_rpe(ii),time(idcs),da(idcs),...
        'color',highlight_clr);
    
    % plot value trace
    stairs(sp_value(ii),time(idcs),value(idcs),...
        'color','k');
    
    % update axis limits
    set(sp_stage(:,ii),...
        'xlim',time_win);
end

% legend
if use_clicks
    leg_str = {'click','reward'};
else
    leg_str = {'reward'};
end
legend(sp_stimulus(1),leg_str,...
    'box','off');

% plot IRI distribution
stem(sp_iri,iri_mu,max([iri_counts,iri_pdf]),...
    'color','k',...
    'marker','v',...
    'markersize',10,...
    'markerfacecolor','k',...
    'markeredgecolor','none',...
    'linewidth',2);
histogram(sp_iri,...
    'binedges',iri_edges,...
    'bincounts',iri_counts .* (iri_edges(1:end-1)>=iri_cutoff),...
    'facecolor','w',...
    'edgecolor','k',...
    'facealpha',1,...
    'linewidth',1);
histogram(sp_iri,...
    'binedges',iri_edges,...
    'bincounts',iri_counts .* (iri_edges(1:end-1)<iri_cutoff),...
    'facecolor',[1,1,1] *.75,...
    'edgecolor','k',...
    'facealpha',1,...
    'linewidth',1);
plot(sp_iri,iri_edges,iri_pdf,...
    'color','k',...
    'linewidth',2);

% plot reaction times
if use_clicks
    plot(sp_reaction,...
        reward_times,reaction_times,...
        'color','k',...
        'marker','.',...
        'markersize',5,...
        'linestyle','none');
end

% plot microstimuli
plot(sp_microstimulus,...
    time(time < max(xlim(sp_microstimulus))),...
    microstimuli(time < max(xlim(sp_microstimulus)),:),...
    'color','k');

% plot eligibility traces
plot(sp_eligibility,...
    time(time < max(xlim(sp_eligibility))),...
    eligibility(time < max(xlim(sp_eligibility)),:),...
    'color','k');

% plot memory trace
plot(sp_microstimulus,...
    time(time < max(xlim(sp_microstimulus))),...
    stimulus_trace(time < max(xlim(sp_microstimulus))),...
    'color','k',...
    'linewidth',2);
text(sp_microstimulus,...
    .05,1,sprintf('y(t) = \\tau^t'),...
    'fontsize',axesopt.fontsize+2,...
    'units','normalized',...
    'horizontalalignment','left');

% highlight a microstimulus
idx = floor(n * 1 / 3);
plot(sp_microstimulus,...
    time(time < max(xlim(sp_microstimulus))),...
    microstimuli(time < max(xlim(sp_microstimulus)),idx),...
    'color','k',...
    'linewidth',2);
plot(sp_eligibility,...
    time(time < max(xlim(sp_eligibility))),...
    eligibility(time < max(xlim(sp_eligibility)),idx),...
    'color','k',...
    'linewidth',2);

% plot average reward-aligned DA signal
for ii = 1 : n_stages
    idcs = (1 : floor(n_rewards/n_stages)) + ...
        (ii -1) * floor(n_rewards/n_stages);
    
    % average reward-aligned reward signal
    da_reward_mu = nanmean(da_reward_snippets(idcs,:));
    plot(sp_reward,da_reward_time,da_reward_mu,...
        'color',stage_clrs(ii,:),...
        'linewidth',1.5);
    
    % average reward-aligned baseline signal
    da_baseline_mu = nanmean(da_baseline_snippets(idcs,:));
    plot(sp_baseline,da_baseline_time,da_baseline_mu,...
        'color',stage_clrs(ii,:),...
        'linewidth',1.5);
end

% legend
legend_labels = repmat({''},n_stages,1);
legend_labels{1} = 'Early';
legend_labels{end} = 'Late';
legend(sp_baseline,legend_labels,...
    'location','northwest',...
    'box','off');

% test 1: DA responses as a function of reward number
nan_flags = isnan(da_response);
reward_idcs = (1 : n_rewards)';
xx = reward_idcs;
[rho,pval] = corr(xx(~nan_flags),da_response(~nan_flags));
mdl = fitlm(xx,da_response);
scatter(sp_test1,xx,da_response,20,...
    colorlerp([early_clr; late_clr],n_rewards),...
    'marker','o',...
    'markerfacealpha',.5,...
    'markerfacecolor','flat',...
    'markeredgecolor','flat',...
    'linewidth',1);
plot(sp_test1,xx,mdl.predict(xx),...
    'color',highlight_clr,...
    'linestyle',repmat('-',1,1+(mdl.Coefficients.pValue(end)>.05)),...
    'linewidth',1.5);
plot(sp_test1,xlim(sp_test1),[1,1]*0,':k');
text(sp_test1,.95,.95,sprintf('r = %.2f',rho),...
    'horizontalalignment','right',...
    'units','normalized',...
    'color','k');
text(sp_test1,.95,.85,sprintf('p = %.2f',pval),...
    'horizontalalignment','right',...
    'units','normalized',...
    'color','k');

% test 2: DA responses as a function of previous IRI
nan_flags = isnan(da_response);
[rho,pval] = corr(iri(~nan_flags),da_response(~nan_flags));
mdl = fitlm(iri,da_response);
scatter(sp_test2,iri,da_response,20,...
    colorlerp([early_clr; late_clr],n_rewards),...
    'marker','o',...
    'markerfacealpha',.5,...
    'markerfacecolor','flat',...
    'markeredgecolor','flat',...
    'linewidth',1);
plot(sp_test2,mdl.Variables.x1,mdl.predict(mdl.Variables.x1),...
    'color',highlight_clr,...
    'linestyle',repmat('-',1,1+(mdl.Coefficients.pValue(end)>.05)),...
    'linewidth',1.5);
plot(sp_test2,xlim(sp_test2),[1,1]*0,':k');
text(sp_test2,.95,.95,sprintf('r = %.2f',rho),...
    'horizontalalignment','right',...
    'units','normalized',...
    'color','k');
text(sp_test2,.95,.85,sprintf('p = %.2f',pval),...
    'horizontalalignment','right',...
    'units','normalized',...
    'color','k');

% axes linkage
arrayfun(@(ax1,ax2,ax3,ax4)linkaxes([ax1,ax2,ax3,ax4],'x'),...
    sp_stimulus,sp_state,sp_rpe,sp_value);
linkaxes([sp_microstimulus,sp_eligibility],'x');
linkaxes(sp_rpe,'y');
linkaxes(sp_value,'y');
linkaxes([sp_baseline,sp_reward],'y');

% annotate model parameters
annotateModelParameters;

%% extra stuff
% [value_snippets,value_time] = ...
%     signal2eventsnippets(time,value,reward_times,[-10,40],dt);
% 
% % figure initialization
% figure(figopt,...
%     'windowstyle','normal',...
%     'name','reward-aligned value');
% set(gca,axesopt,...
%     'plotboxaspectratio',[2,1,1]);
% 
% % axes labels
% xlabel('Time since reward (s)');
% ylabel('Value (a.u.)');
% 
% % graphical object preallocation
% p = gobjects(n_stages,1);
% 
% % iterate through training stages
% for ii = 1 : n_stages
%     idcs = (1 : floor(n_rewards/n_stages)) + ...
%         (ii -1) * floor(n_rewards/n_stages);
%     errorpatch(...
%         value_time,...
%         nanmean(value_snippets(idcs,:)),...
%         nanstd(value_snippets(idcs,:))./sqrt(numel(idcs)),...
%         stage_clrs(ii,:),...
%         'facealpha',.25);
%     p(ii) = plot(value_time,nanmean(value_snippets(idcs,:)),...
%         'linewidth',1.5,...
%         'color',stage_clrs(ii,:));
% end
% 
% % alignemnt line
% plot([0,0],ylim,'--k');
% 
% % legend
% legend_labels = repmat({''},n_stages,1);
% legend_labels{1} = 'Early';
% legend_labels{end} = 'Late';
% legend(p,legend_labels,...
%     'location','northeast',...
%     'box','off');