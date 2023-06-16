%% initialization
close all;
clear;
clc;

%% directory settings
root_path = fileparts(pwd);
data_path = fullfile(root_path,'data');
data_dir = dir(data_path);
data_dir = data_dir(cellfun(@(x)~contains(x,'.'),{data_dir.name}));

%% mouse settings
mouse_ids = {data_dir.name};
n_mice = numel(mouse_ids);

%% experiment settings
experiment_id = 'RandomRewards';

%% analysis parameters
baseline_period = [-2,-.5];
event_period = [-.5,1];
iri_cutoff = 3;

%% data parsing

% iterate through mice
for mm = 1 : n_mice
    progressreport(mm,n_mice,'parsing behavioral data');
    
    mouse_path = fullfile(data_path,mouse_ids{mm},experiment_id);
    mouse_dir = dir(mouse_path);
    mouse_dir = mouse_dir(cellfun(@(x)~contains(x,'.'),{mouse_dir.name}));
    
    session_ids = {mouse_dir.name};
    session_days = cellfun(@(x)str2double(strrep(x,'Day','')),session_ids);
    [~,sorted_idcs] = sort(session_days);
    session_ids = session_ids(sorted_idcs);
    n_sessions = numel(session_ids);
    
    % initialize mouse counters
    mouse_dur_counter = 0;
    mouse_trial_counter = 0;
    mouse_lick_counter = 0;
    
    % iterate through sessions
    for ss = 1 : n_sessions
        session_path = fullfile(mouse_path,session_ids{ss});
        
        %% behavior
        
        % load behavior
        bhv_id = sprintf('%s_%s_eventlog.mat',...
            mouse_ids{mm},session_ids{ss});
        bhv_path = fullfile(session_path,bhv_id);
        load(bhv_path)
        
        % parse behavior
        event_labels = categorical(...
            eventlog(:,1),[5,7,0],{'lick','reward','end'});
        event_times = eventlog(:,2);
        event_idcs = (1 : numel(event_labels))';
        start_idx = find(event_labels == 'reward',1);
        end_idx = find(event_labels == 'end');
        session_dur = event_times(end_idx);
        
        valid_flags = ...
            ~isundefined(event_labels) & ...
            event_idcs >= start_idx & ...
            event_idcs < end_idx;
        n_events = sum(valid_flags);
        event_idcs = (1 : n_events)';
        event_labels = removecats(event_labels(valid_flags),'end');
        event_times = event_times(valid_flags);
        
        reward_flags = event_labels == 'reward';
        reward_idcs = find(reward_flags);
        reward_times = event_times(reward_flags);
        n_rewards = sum(reward_flags);
        
        lick_flags = event_labels == 'lick';
        lick_trials = nan(n_events,1);
        
        % iterate through rewards
        for ii = 1 : n_rewards
            trial_onset = reward_idcs(ii);
            if ii < n_rewards
                trial_offset = reward_idcs(ii+1);
            else
                trial_offset = inf;
            end
            trial_flags = ...
                event_idcs > trial_onset & ...
                event_idcs < trial_offset;
            n_licks = sum(lick_flags & trial_flags);
            lick_trials(trial_flags) = 1 : n_licks;
        end
        
        lick_sessions = cumsum(~isnan(lick_trials),'omitnan');
        lick_sessions(lick_sessions == 0) = nan;
        
        solenoid_times = event_times(reward_flags);
        event_trials = sum(event_times >= event_times(reward_flags)',2);
        n_trials = max(event_trials);
        iri_prev = diff([0;solenoid_times]);
        iri_next = diff([solenoid_times;nan]);
        
        %% photometry
        photometry_path = fullfile(session_path,'Photometry.mat');
        load(photometry_path);
        
        % correct for nonstationary sampling frequency
        fs = 120;
        dt = 1 / fs;
        time = T(1) : dt : T(end);
        da = interp1(T,dff,time);
        
        % get event-aligned snippets of DA
        da_baseline_snippets = ...
            signal2eventsnippets(time,da,event_times,baseline_period,dt);
        da_event_snippets = ...
            signal2eventsnippets(time,da,event_times,event_period,dt);
        
        % preallocation
        da_response = nan(n_events,1);
        
        % iterate through events
        for ii = 1 : n_events
            
            % compute 'DA response' metric
            da_response(ii) = ...
                sum(da_event_snippets(ii,:)) - sum(da_baseline_snippets(ii,:));
        end
        
        %% combine behavioral & photometry data
        trial_subtable = table(...
            event_trials + mouse_trial_counter,...
            event_trials,...
            'variablenames',{'mouse','session'});
        
        lick_subtable = table(...
            lick_sessions + mouse_lick_counter,...
            lick_sessions,...
            lick_trials,...
            'variablenames',{'mouse','session','trial'});
        
        time_subtable = table(...
            event_times+mouse_dur_counter,...
            event_times,...
            event_times-solenoid_times(event_trials),...
            'variablenames',{'mouse','session','trial'});
        
        iri_table = table(...
            iri_prev(event_trials),...
            iri_next(event_trials),...
            'variablenames',{'previous','next'});
        
        da_table = table(...
            da_baseline_snippets,...
            da_event_snippets,...
            da_response,...
            'variablenames',{'pre','post','response'});
        
        session_events = table(...
            categorical(cellstr(repmat(mouse_ids{mm},n_events,1)),mouse_ids),...
            repmat(ss,n_events,1),...
            event_labels,...
            trial_subtable,...
            lick_subtable,...
            time_subtable,...
            iri_table,...
            da_table,...
            'variablenames',...
            {'mouse','session','label','trial','lick','time','iri','DA'});
        
        %         figure('windowstyle','docked');
        %
        %         subplot(2,1,1);
        %         plot(event_times-reward_times(event_trials),event_trials,'.k');
        %         title(bhv_id,'interpreter','none');
        %         xlabel('Time since reward (s)');
        %         xlim([0,12])
        %
        %         subplot(2,1,2);
        %         reaction_times = event_times(reward_idcs+1) - reward_times;
        %         histogram(reaction_times,linspace(0,5,50),...
        %             'facecolor','k');
        %         xlabel('Reaction time (s)');
        
        % increment mouse counter
        mouse_dur_counter = mouse_dur_counter + session_dur;
        mouse_trial_counter = mouse_trial_counter + n_trials;
        mouse_lick_counter = mouse_lick_counter + ...
            lick_sessions(find(~isnan(lick_sessions),1,'last'));
        
        if ss == 1
            mouse_events = session_events;
        else
            mouse_events = [mouse_events;session_events];
        end
    end
    
    if mm == 1
        events = mouse_events;
    else
        events = [events;mouse_events];
    end
end

%%

% iterate through mice
for mm = 1 : n_mice
    mouse_flags = events.mouse == mouse_ids{mm};
    lick_flags = mouse_flags & events.label == 'lick';
    reward_flags = mouse_flags & lick_flags & events.lick.trial == 1;
    
    lick_times = events.time.mouse(lick_flags);
    ili = [nan;diff(lick_times)];
    
    reward_idcs = find(reward_flags);
    
    reaction_times = ...
        events.time.session(reward_idcs) - ...
        events.time.session(reward_idcs-1);
    
    trial_iri = events.iri.previous(reward_idcs-1);
    iri_flags = ...
        [trial_iri(1);trial_iri(1:end-1)] >= 3;
    
    
    figure('windowstyle','docked');
    
    subplot(4,2,[1,3]);
    hold on;
    plot(events.time.trial(lick_flags),...
        events.trial.mouse(lick_flags),'.k',...
        'markersize',5);
    plot(events.time.trial(reward_idcs),...
        events.trial.mouse(reward_idcs),'or');
    title(mouse_ids{mm},'interpreter','none');
    xlabel('Time since reward (s)');
    ylabel('Trial #');
    axis tight;
    xlim([0,3]);
    
    subplot(4,2,5);
    hold on;
    set(gca,...
        'ytick',1:100:1e3);
    histogram(reaction_times,linspace(0,5,100),...
        'facecolor','k');
    histogram(reaction_times(iri_flags),linspace(0,5,100),...
        'facecolor','r');
    xlabel('Reaction time (s)');
    axis tight;
    
    subplot(4,2,7);
    plot(events.trial.mouse(reward_idcs),reaction_times,'.k');
    axis tight;
    ylim([0,quantile(reaction_times,.975)]);
    xlabel('Trial #');
    ylabel('Reaction time (s)');
    
    subplot(4,2,[2,4]);
    plot(events.trial.mouse(lick_flags),ili,'.k',...
        'markersize',5);
    xlabel('Trial #');
    ylabel('Inter-lick-interval (s)');
    axis tight;
    ylim(quantile(ili,[.01,.8]));
    
    subplot(4,2,[6,8]);
    plot(events.trial.mouse(lick_flags),...
        events.iri.previous(lick_flags),'.k',...
        'markersize',5);
    xlabel('Trial #');
    ylabel('Inter-reward-interval (s)');
    axis tight;
    
    valid_flags = ...
        events.iri.previous > iri_cutoff & ...
        events.iri.next > iri_cutoff;

    % replot fig3E
    figure('color','w'); hold on;
    set(gca,...
        'linewidth',2,...
        'fontsize',12,...
        'tickdir','out');
    title(mouse_ids{mm},...
        'interpreter','none');
    xlabel('Reward #');
    ylabel('DA response');
    n_sessions = max(events.session(mouse_flags));
    clrs = cool(n_sessions);
    for ss = 1 : n_sessions
        session_flags = events.session == ss;
        trial_flags = ...
            mouse_flags & ...
            session_flags & ...
            reward_flags;
        x = events.trial.mouse(trial_flags & valid_flags);
        X = [ones(size(x)),x];
        y = events.DA.response(trial_flags & valid_flags);
        betas = robustfit(x,y);
        plot(events.trial.mouse(trial_flags),events.DA.response(trial_flags),'.',...
            'markersize',10,...
            'color',[1,1,1]*.85);
        plot(x,y,'.',...
            'markersize',15,...
            'color',clrs(ss,:));
        plot(x,X*betas,'-w',...
            'linewidth',3);
        plot(x,X*betas,'-',...
            'color',clrs(ss,:),...
            'linewidth',1.5);
    end
    trial_flags = ...
        mouse_flags & ...
        reward_flags;
    x = events.trial.mouse(trial_flags & valid_flags);
    X = [ones(size(x)),x];
    y = events.DA.response(trial_flags & valid_flags);
    nan_flags = isnan(y);
    betas = robustfit(x,y);
    plot(x,X*betas,'--k',...
        'linewidth',1);
    axis tight;
end