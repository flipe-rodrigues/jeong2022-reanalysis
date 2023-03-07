function [snippets,alignment_time] = signal2eventsnippets(...
    time,signal,event_times,alignment_window,delta_time)
    %SIGNAL2EVENTSNIPPETS Takes a SIGNAL as input and outputs a matrix of 
    %signal SNIPPETS spanning ALIGNMENT_WINDOW around EVENT_TIMES.

    % input parsing
    n_events = numel(event_times);
    alignment_window = delta_time * round(alignment_window / delta_time);
    alignment_time = alignment_window(1) : delta_time : alignment_window(2);
    n_samples = numel(alignment_time);
    duration = range(time) + delta_time;

    % preallocation
    snippets = nan(n_events,n_samples);

    % iterate through events
    for ii = 1 : n_events
        progressreport(ii,n_events,'aligning');
        onset_time = event_times(ii) + alignment_window(1);
        offset_time = event_times(ii) + alignment_window(2);
        snippet_time = onset_time : delta_time : offset_time;
        snippet_flags = snippet_time >= 0 & snippet_time < duration;
        [~,onset_idx] = min(abs(time - onset_time));
        [~,offset_idx] = min(abs(time - offset_time));
        idcs = onset_idx : offset_idx;
        snippets(ii,snippet_flags) = signal(idcs);
    end
end

