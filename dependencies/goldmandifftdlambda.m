function [state,value,rpe,rwdrate,weights,eligibility] = goldmandifftdlambda(...
        time,...
        stimulus_presence,...
        reward_times,...
        weights,...
        varargin)
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here

    %% input parsing
    n_states = numel(time);
    n_stimuli = size(stimulus_presence,2);
    dt = diff(time(1:2));
    duration = n_states * dt;

    % name-value pair input parsing
    p = inputParser;
    p.addParameter('gamma',.98);
    p.addParameter('alpha',.01);
    p.addParameter('lambda',.95);
    p.addParameter('tau',.95);
    p.addParameter('theta',0);
    p.addParameter('n',20);
    p.parse(varargin{:});
    param = p.Results;

    %% network parameters
    tau = .25; % in seconds
    Wi = normpdf(1:param.n,(1:param.n)',param.n);
    Wi = -(1 - normalize01(Wi)) * 1 / param.n * 0;
    W = circshift(eye(param.n),-1) .* [ones(param.n-1,1);0] + Wi;
    I = eye(1,param.n) - (ones(1,param.n) - eye(1,param.n)) * .075;
    
    %% construct state features
    
    % preallocation
    stimulus_features = zeros(n_states,param.n,n_stimuli);

    % iterate through stimuli
    for ii = 1 : n_stimuli

        % input corresponding to the current stimulus
        x = stimulus_presence(:,ii);
        
        % preallocation
        r = zeros(n_states,param.n);
        
        % iterate through states
        for ss = 2 : n_states
            progressreport(ss-1,n_states-1,sprintf(...
                'constructing stimulus features (%i/%i)',ii,n_stimuli));
            drdt = (-r(ss-1,:) + r(ss-1,:) * W + x(ss) * I) / (tau / dt);
            r(ss,:) = max(0,r(ss-1,:) + drdt);
        end
        
        % store stimulus features
        stimulus_features(:,:,ii) = r * 2.5; % r / max(r(:)) * 2;
    end

    % concatenate representations across stimuli
    state = reshape(stimulus_features,n_states,n_stimuli*param.n);
    n_features = size(state,2);

    %% compute reward function
    time_edges = linspace(0,duration,n_states+1);
    reward = histcounts(reward_times,time_edges);
%     reward = state(:,1);
    
    %% differential td learning
    
    % preallocation
    value = zeros(n_states,1);
    rpe = zeros(n_states,1);
    rwdrate = zeros(n_states,1);
    eligibility = zeros(n_states,n_features);
    if isempty(weights)
        weights = zeros(1,n_features);
    end
    
    % iterate through states
    for ss = 2 : n_states
        progressreport(ss-1,n_states-1,'running differential TD(lambda)');
        eligibility(ss,:) = ...
            param.lambda * eligibility(ss-1,:) + state(ss-1,:);
        value(ss) = weights * state(ss,:)';
        rpe(ss) = reward(ss) - rwdrate(ss - 1) + value(ss) - value(ss-1);
%         rpe(ss) = reward(ss) + param.gamma * value(ss) - value(ss-1);
        weights = weights + param.alpha * rpe(ss) * eligibility(ss,:) * dt;
        rwdrate(ss) = rwdrate(ss-1) + ...
            (1 - param.gamma) * 10 * param.alpha * rpe(ss);
    end
    
    % convert to units of time
    rwdrate = rwdrate / dt;
end