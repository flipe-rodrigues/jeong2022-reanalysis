function [snippets,alignment_time] = signal2eventsnippets(...
    time,signal,event_times,alignment_window,delta_time,nanify)
    %SIGNAL2EVENTSNIPPETS Takes a SIGNAL as input and outputs a matrix of
    %signal SNIPPETS spanning ALIGNMENT_WINDOW around EVENT_TIMES.

    if nargin < 6
        nanify = false;
    end
    
    % input parsing
    n_events = numel(event_times);
    alignment_window = delta_time * round(alignment_window / delta_time);
    alignment_time = alignment_window(1) : delta_time : alignment_window(2);
    n_samples = numel(alignment_time);

    % preallocation
    snippets = nan(n_events,n_samples);

    % iterate through events
    for ii = 1 : n_events
        onset_time = event_times(ii) + alignment_window(1);
        offset_time = event_times(ii) + alignment_window(2);
        snippet_time = onset_time : delta_time : offset_time;
        snippet_flags = ...
            snippet_time > time(1) - delta_time / 2 & ...
            snippet_time < time(end) + delta_time / 2;
        time_flags = ...
            time > onset_time - delta_time / 2 & ...
            time < offset_time + delta_time / 2;
        snippets(ii,snippet_flags) = signal(time_flags);
    end
    
    % nanify overlapping snippets
    if nanify
        iei = diff([0;event_times]);
        valid_mask = ...
            alignment_time > -iei & ...
            alignment_time < [iei(2:end);inf];
        snippets(~valid_mask) = nan;
    end
end