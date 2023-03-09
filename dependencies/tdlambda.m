function [state,value,rpe,weights,eligibility] = tdlambda(...
        time,stimulus_times,reward_times,temporal_bases,weights,varargin)
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here

    %% input parsing
    n_stimuli = size(stimulus_times,2);
    [n_states,n_bases] = size(temporal_bases);

    % name-value pair input parsing
    p = inputParser;
    p.addParameter('gamma',.98);
    p.addParameter('alpha',.01);
    p.addParameter('lambda',.95);
    p.addParameter('tau',.95);
    p.addParameter('theta',0);
    p.parse(varargin{:});
    param = p.Results;

    %% construct state features
    
    % preallocation
    stimulus_features = zeros(n_states,n_bases,n_stimuli);

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
                'constructing state features');
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
            stimulus_features(idcs,:,ii) = temporal_bases(1:horizon,:);
        end
    end

    % concatenate representations across stimuli
    state = reshape(stimulus_features,n_states,n_stimuli*n_bases);

    %% compute reward function
    dt = diff(time(1:2));
    duration = n_states * dt;
    time_edges = linspace(0,duration,n_states+1);
    reward = histcounts(reward_times,time_edges);
    
    %% td learning
    
    % preallocation
    value = zeros(n_states,1);
    rpe = zeros(n_states,1);
    eligibility = zeros(n_states,n_stimuli*n_bases);
    z = zeros(1,n_stimuli*n_bases);
    if isempty(weights)
        weights = zeros(1,n_stimuli*n_bases);
    end
    
    % iterate through states
    for ss = 2 : n_states
        progressreport(ss-1,n_states-1,'running TD(lambda)');
        z = param.gamma * param.lambda * z + state(ss-1,:);
        eligibility(ss,:) = z;
        value(ss) = weights * state(ss,:)';
        rpe(ss) = reward(ss) + param.gamma * value(ss) - value(ss-1);
        weights = weights + param.alpha * rpe(ss) * z * dt;
    end
end