% iterate through mice
for mm = 1 : n_mice
    mouse_flags = data.mouse == mouse_ids{mm};
    reward_flags = ...
        valid_flags & ...
        mouse_flags;
    n_rewards = sum(reward_flags);
    if n_rewards == 0
        continue;
    end
    
    % parse session data
    n_sessions = max(data.session(mouse_flags));
    n_rewards_per_session = arrayfun(@(x) ...
        sum(data.session(reward_flags) == x),1:n_sessions);
    session_delims = [1,cumsum(n_rewards_per_session)];
    session_clrs = cool(n_sessions);
    
    % figure initialization
    figure(...
        'numbertitle','off',...
        'name',sprintf('%s_%s_da_rasters_%s_%s',...
        meta.experiment,mouse_ids{mm},alignment,strjoin(sorting_criteria,'_')),...
        'color','w');
    
    % axes initialization
    yytick = session_delims;
    yyticklabel = num2cell(yytick);
    yyticklabel(2:end-1) = {''};
    axes(axesopt.default,...
        'xlim',[-1,1]*5,...
        'ylim',[0,n_rewards],...
        'ytick',yytick,...
        'yticklabel',yyticklabel,...
        'colormap',bone(2^8-1));
    
    % axes labels
    title(sprintf('%s',mouse_ids{mm}),...
        'interpreter','none');
    xlabel(sprintf('Time since reward %s (s)',alignment));
    ylabel(sprintf('Reward # (sorted by %s)',sorting_criteria{end}));
    
    % reward sorting
    [~,sorted_idcs] = sortrows(data(reward_flags,:),sorting_criteria);
    
    % plot DA raster
    da_mat = data.da.(alignment)(reward_flags,:);
    imagesc(meta.epochs.(alignment),[0,n_rewards]+[1,-1]*.5,...
        da_mat(sorted_idcs,:),quantile(da_mat,[.001,.999],'all')');
    
    % plot reaction times
    reaction_times = data.rt(reward_flags);
    session_idcs = data.session(reward_flags);
    scatter(...
        reaction_times(sorted_idcs),1:n_rewards,7.5,...
        session_clrs(session_idcs(sorted_idcs),:),...
        'marker','.');
    
    % iterate through sessions
    for ss = 1 : n_sessions
        session_flags = data.session == ss;
        reward_flags = ...
            valid_flags & ...
            mouse_flags & ...
            session_flags;
        
        % increment reward counters
        prev_counter = curr_counter;
        curr_counter = curr_counter + sum(reward_flags);
        
        % plot session delimeters
        if ss < n_sessions
            plot(xlim,[1,1]*session_delims(ss+1),...
                'color',[1,1,1],...
                'linestyle',':');
        end
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
    
    % plot reference lines
    plot([0,0],ylim,'--w');
    
    % save figure
    if want2savepanels
        file_name = sprintf('%s',[get(gcf,'name'),'.png']);
        file_path = fullfile(panel_path,file_name);
        print(gcf,file_path,'-dpng','-r300','-painters');
    end
end