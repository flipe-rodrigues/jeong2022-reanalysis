function [state,value,rpe,weights,eligibility] = inhibitiontdlambda(...
        time,...
        presence_features,...
        stimulus_times,...
        reward_times,...
        temporal_bases,...
        weights,...
        rpe_inhibition,...
        varargin)
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here

    %% input parsing
    n_states = numel(time);
    n_stimuli = size(stimulus_times,2);
    n_bases = size(temporal_bases,2);

    % name-value pair input parsing
    p = inputParser;
    p.addParameter('gamma',.98);
    p.addParameter('alpha',.01);
    p.addParameter('lambda',.95);
    p.addParameter('tau',.95);
    p.parse(varargin{:});
    param = p.Results;

    %% construct state features
    
    % preallocation
    stimulus_temporal_features = zeros(n_states,n_bases,n_stimuli);

    % for progress keeping purposes
    iteration_count = sum(cellfun(@(x)sum(~isnan(x)),stimulus_times));
    iteration = 1;
    
    % iterate through stimuli
    for ii = 1 : n_stimuli
        
        % get current stimulus times
        current_stimulus_times = stimulus_times{:,ii};
        nan_flags = isnan(current_stimulus_times);
        [~,stimulus_state_idcs] = ...
            min(abs(time - current_stimulus_times(~nan_flags)),[],2);
        stimulus_count = numel(stimulus_state_idcs);
        
        % iterate through stimuli
        for jj = 1 : stimulus_count
            progressreport(iteration,iteration_count,...
                'constructing stimulus features');
            iteration = iteration + 1;
            if jj == stimulus_count
                horizon = n_states - stimulus_state_idcs(jj);
            else
                horizon = stimulus_state_idcs(jj+1) - stimulus_state_idcs(jj);
            end
            idcs = (1 : horizon) + stimulus_state_idcs(jj) - 1;
%             temporal_scaling = max(0,normrnd(1,param.theta));
%             unscaled_bases = temporal_bases;
%             if horizon > 1
%                 scaled_bases = interp1(...
%                     time,unscaled_bases,time*temporal_scaling);
%             else
%                 scaled_bases = unscaled_bases;
%             end
            stimulus_temporal_features(idcs,:,ii) = ...
                ...1/2 * stimulus_temporal_features(idcs,:,ii) + ...
                temporal_bases(1:horizon,:);
        end
    end

    % concatenate representations across stimuli
    state = reshape(stimulus_temporal_features,n_states,n_stimuli*n_bases);
    
    %% concatenate stimulus & temporal features
    presence_features = ...
        normalize01(presence_features,1) * max(temporal_bases,[],'all');
    state = [state,presence_features];
    n_features = size(state,2);
    
    %% compute reward function
    dt = diff(time(1:2));
    duration = n_states * dt;
    time_edges = linspace(0,duration,n_states+1);
    reward = histcounts(reward_times,time_edges);
    
    %% td learning
    
    % preallocation
    value = zeros(n_states,1);
    rpe = zeros(n_states,1);
    eligibility = zeros(n_states,n_features);
    if isempty(weights)
        weights = zeros(1,n_features);
    end
    
    % iterate through states
    for ss = 2 : n_states
        progressreport(ss-1,n_states-1,'running TD(lambda)');
        eligibility(ss,:) = ...
            param.gamma * param.lambda * eligibility(ss-1,:) + state(ss-1,:);
        value(ss) = weights * state(ss,:)';
%         if rpe_inhibition(ss) == 0
%             rpe(ss) = reward(ss) + param.gamma * value(ss) - value(ss-1);
%         else
%             rpe(ss) = rpe_inhibition(ss);
%         end
        rpe(ss) = reward(ss) + param.gamma * value(ss) - value(ss-1) + ...
            rpe_inhibition(ss);
        weights = weights + param.alpha * rpe(ss) * eligibility(ss,:) * dt;
    end
end