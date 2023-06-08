%% initialization
close all;
clear;
clc;

%%
root_path = fileparts(pwd);
data_path = fullfile(root_path,'data');
data_dir = dir(data_path);
data_dir = data_dir(cellfun(@(x)~contains(x,'.'),{data_dir.name}));

%%
mouse_ids = {data_dir.name};
n_mice = numel(mouse_ids);

%%
experiment_id = 'RandomRewards';

%%

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
    
    % initialize mouse trial counter
    mouse_trial_counter = 0;
    
    % iterate through sessions
    for ss = 1 : n_sessions
        
        session_path = fullfile(mouse_path,session_ids{ss});
        bhv_id = sprintf('%s_%s_eventlog.mat',...
            mouse_ids{mm},session_ids{ss});
        bhv_path = fullfile(session_path,bhv_id);
        load(bhv_path)
        
        event_labels = categorical(...
            eventlog(:,1),[5,7,0],{'lick','reward','end'});
        event_times = eventlog(:,2);
        event_idcs = (1 : numel(event_labels))';
        start_idx = find(event_labels == 'reward',1);
        end_idx = find(event_labels == 'end');
        
        valid_flags = ...
            ~isundefined(event_labels) & ...
            event_idcs >= start_idx & ...
            event_idcs < end_idx;
        n_events = sum(valid_flags);
        event_labels = removecats(event_labels(valid_flags),'end');
        event_times = event_times(valid_flags);
        reward_flags = event_labels == 'reward';
        reward_times = event_times(reward_flags);
        event_trials = sum(event_times >= event_times(reward_flags)',2);
        n_trials = max(event_trials);
        iri = diff([0;reward_times]);
        
        trial_subtable = table(...
            event_trials + mouse_trial_counter,...
            event_trials,...
            'variablenames',{'mouse','session'});
        
        time_subtable = table(...
            event_times,...
            event_times-reward_times(event_trials),...
            'variablenames',{'session','trial'});
        
        session_events = table(...
            categorical(cellstr(repmat(mouse_ids{mm},n_events,1)),mouse_ids),...
            repmat(ss,n_events,1),...
            event_labels,...
            trial_subtable,...
            time_subtable,...
            iri(event_trials),...
            'variablenames',{'mouse','session','label','trial','time','iri'});
        
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
        
        % increment mouse trial counter
        mouse_trial_counter = mouse_trial_counter + n_trials;
        
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

    figure('windowstyle','docked');
    
    subplot(4,1,[1,2]);
    lick_flags = mouse_events.label == 'lick';
    plot(mouse_events.time.trial(lick_flags),...
        mouse_events.trial.mouse(lick_flags),'.k');
    title(mouse_ids{mm},'interpreter','none');
    xlabel('Time since reward (s)');
    axis tight;
    xlim([0,12]);
    
    subplot(4,1,3);
    hold on;
    set(gca,...
        'ytick',1:100:1e3);
    reward_flags = mouse_events.label == 'reward';
    reward_idcs = find(reward_flags);
    reward_times = mouse_events.time.session(reward_idcs);
    lick_times = mouse_events.time.session(mouse_events.label == 'lick');
    reaction_times = ...
        mouse_events.time.session(reward_idcs+1) - ...
        mouse_events.time.session(reward_idcs);
    trial_iri = mouse_events.iri(reward_idcs);
    iri_flags = ...
        [trial_iri(1);trial_iri(1:end-1)] >= 3;
    histogram(reaction_times,linspace(0,5,100),...
        'facecolor','k');
    histogram(reaction_times(iri_flags),linspace(0,5,100),...
        'facecolor','r');
    xlabel('Reaction time (s)');
    axis tight;
    
    subplot(4,1,4);
    plot(mouse_events.trial.mouse(reward_idcs),reaction_times,'.k');
    axis tight;
    ylim([0,quantile(reaction_times,.975)]);
end