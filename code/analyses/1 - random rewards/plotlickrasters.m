function plotlickrasters(data,meta,opt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% argument validation
arguments
    data table
    meta struct
    opt.alignment string
    opt.sorting (:,:) cell
    opt.xlim (1,2) double
    opt.xlabel string
    opt.ylabel string
    opt.savepath string
end

% iterate through mice
for mm = 1 : meta.mice.n
    mouse_flags = data.mouse == meta.mice.ids{mm};
    trial_flags = ...
        mouse_flags;
    n_trials = sum(trial_flags);
    
    % parse session data
    n_sessions = max(data.session(mouse_flags));
    n_rewards_per_session = arrayfun(@(x) ...
        sum(data.session(trial_flags) == x),1:n_sessions);
    session_delims = [1,cumsum(n_rewards_per_session)];
    session_clrs = cool(n_sessions);
    
    % figure initialization
    figure(...
        'numbertitle','off',...
        'name',sprintf('%s_%s_lick_rasters_%s_%s',...
        meta.experiment,meta.mice.ids{mm},opt.alignment,opt.sorting{end}),...
        'color','w');
    
    % axes initialization
    yytick = session_delims;
    yyticklabel = num2cell(yytick);
%     yyticklabel(2:end-1) = {''};
    axes(defaultaxesoptions,...
        'xlim',opt.xlim,...
        'ylim',[0,n_trials],...
        'ytick',yytick,...
        'yticklabel',yyticklabel,...
        'colormap',bone(2^8-1));
    
    % axes labels
    title(sprintf('%s',meta.mice.ids{mm}));
    xlabel(opt.xlabel);
    ylabel(opt.ylabel);
    
    % trial sorting
    [~,sorted_idcs] = sortrows(data(trial_flags,:),opt.sorting);
    
    % fetch lick times
    lick_times = data.lick.time.(opt.alignment)(trial_flags);
    
    % iterate through trials
    for tt = 1 : n_trials
        plot(lick_times{sorted_idcs(tt)},...
            repmat(tt,numel(lick_times{sorted_idcs(tt)}),1),...
            'color','k',...
            'marker','.',...
            'markersize',3,...
            'linestyle','none');
    end
    
    % plot reaction times
    rt_sign = (-1) ^ strcmpi(opt.alignment,'collection');
    rt = data.rt(trial_flags) * rt_sign;
    session_idcs = data.session(trial_flags);
    scatter(...
        rt(sorted_idcs),1:n_trials,10,...
        session_clrs(session_idcs(sorted_idcs),:),...
        'marker','.');
    
    % iterate through sessions
    for ss = 1 : n_sessions
        
        % plot session delimeters
        if ss < n_sessions
%             plot(xlim,[1,1]*session_delims(ss+1),...
%                 'color','k',...
%                 'linestyle',':');
        end
        
        % patch session edge bands
        xpatch = ...
            [1,1,1,1] * min(xlim) + ...
            [0,1,1,0] * range(xlim) * .025;
        ypatch = ...
            [1,1,0,0] * session_delims(ss) + ...
            [0,0,1,1] * session_delims(ss+1);
        patch(xpatch,ypatch,session_clrs(ss,:),...
            'edgecolor','none',...
            'facealpha',1);
    end
    
    % plot alignment line
%     plot([0,0],ylim,'--k');
    
    % save figure
    if ~isempty(opt.savepath)
        file_name = sprintf('%s',[get(gcf,'name'),'.png']);
        file_path = fullfile(opt.savepath,file_name);
        print(gcf,file_path,'-dpng','-r300','-painters');
    end
end
end