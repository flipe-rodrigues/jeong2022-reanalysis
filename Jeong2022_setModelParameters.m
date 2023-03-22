%% converts from/to temporal discount factor to/from time constant [s]
gammafun = @(dt,eta) exp(-dt/eta);
etafun = @(dt,gamma) -dt/log(gamma);

%% model parameters (partially inspired by Wei et al. 2022)
dt = .02;                               % state step size [s]
eta = 10;                            	% discount time constant [s]
gamma = gammafun(dt,eta);               % temporal discount factor
alpha = .15;                            % learning rate
lambda = gammafun(.02,etafun(.2,.95));	% decay for eligibility traces
tau = .95;                              % decay for the stimulus trace
theta = .15;                            % std of temporal scaling (NOT IN USE!!!)
y0 = 1;                                 % starting height of the stimulus trace
sigma = .08;                            % width of each basis function
n = 20;                                 % number of microstimuli per stimulus

%% model parameters (high resolution version of Jeong & Namboodiri 2022)
% dt = .02;                               % state step size [s]
% eta = etafun(.2,.98);                   % discount time constant [s]
% gamma = gammafun(dt,eta);               % temporal discount factor
% alpha = .02;                            % learning rate
% lambda = gammafun(.02,etafun(.2,.95));	% decay for eligibility traces
% tau = .99;                              % decay for the stimulus trace
% y0 = 1;                                 % starting height of the stimulus trace
% sigma = .08;                            % width of each basis function
% n = 20;                                 % number of microstimuli per stimulus

%% model parameters (matched to Ludvig & Sutton 2008)
% dt = .05;                   % state step size [s]
% gamma = .98;                % temporal discount factor
% eta = etafun(dt,gamma);     % discount time constant [s]
% alpha = .01;                % learning rate
% lambda = .95;               % decay for eligibility traces
% tau = .985;                 % decay for the stimulus trace
% y0 = 1;                     % starting height of the stimulus trace
% sigma = .08;                % width of each basis function
% n = 50;                     % number of microstimuli per stimulus

%% model parameters (matched to Jeong & Namboodiri 2022)
% dt = .2;                    % state step size [s]
% gamma = .98;                % temporal discount factor
% eta = etafun(dt,gamma);     % discount time constant [s]
% alpha = .02;                % learning rate
% lambda = .95;               % decay for eligibility traces
% tau = .99;                  % decay for the stimulus trace
% y0 = 1;                     % starting height of the stimulus trace
% sigma = .08;                % width of each basis function
% n = 20;                     % number of microstimuli per stimulus