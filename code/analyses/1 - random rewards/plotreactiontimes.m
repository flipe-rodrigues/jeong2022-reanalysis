function plotreactiontimes(data,meta,opt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% argument validation
arguments
    data table
    meta struct
    opt.alignment string
    opt.xlim (1,2) double
    opt.xlabel string
    opt.savepath string
end

%
entropyfun = @(p) -nansum(p .* log2(p));

% reaction time bin settings
rt_binwidth = 1 / 15;
rt_roi = meta.epochs.(opt.alignment);
rt_edges = rt_roi(1) : rt_binwidth : rt_roi(2);

% iterate through mice
for mm = 1 : meta.mice.n
    mouse_flags = data.mouse == meta.mice.ids{mm};
    
    % parse session data
    n_sessions = max(data.session(mouse_flags));
    session_clrs = cool(n_sessions);
    
    % figure initialization
    figure(...
        'numbertitle','off',...
        'name',sprintf('%s_%s_rt_%s',...
        meta.experiment,meta.mice.ids{mm},opt.alignment),...
        'color','w');
    
    % axes initialization
    axes(defaultaxesoptions,...
        'plotboxaspectratio',[1,1,1],...
        'xlim',opt.xlim,...
        'ycolor','none');
    
    % axes labels
    title(sprintf('%s',meta.mice.ids{mm}));
    xlabel(opt.xlabel);
    
    % iterate through sessions
    for ss = 1 : n_sessions
        session_flags = data.session == ss;
        trial_flags = ...
            mouse_flags & ...
            session_flags;
        
        % vertical offset
        offset = -ss;
        
        % compute reaction time distribution
        rt = data.rt(trial_flags);
        rt_counts = histcounts(rt,rt_edges);
        rt_x = rt_edges(1:end-1);
        rt_y = normalize01(rt_counts,2) * .85;
        rt_x_ups = sort([rt_edges(2:end),rt_edges(2:end)]);
        rt_y_ups = upsample(rt_y,2);
        rt_y_ups(2:2:end-1) = rt_y(2:end);
        xpatch = [rt_x_ups,fliplr(rt_x_ups)];
        ypatch = [rt_y_ups,zeros(size(rt_y_ups))];
        patch(xpatch,ypatch+offset,session_clrs(ss,:),...
            'edgecolor','none',...
            'facecolor',session_clrs(ss,:),...
            'facealpha',.25,...
            'linewidth',1.5);
        stairs(rt_x,rt_y+offset,...
            'color',session_clrs(ss,:),...
            'linewidth',1.5);
        
        % compute basal entropy
        iri_pdf = histcounts(data.iri.nominal(trial_flags),...
            'binwidth',rt_binwidth,...
            'normalization','pdf');
        iri_pdf(iri_pdf == 0) = nan;
        hb = entropyfun(iri_pdf);
        
        % compute residual entropy
        rt_pdf = rt_counts / nansum(rt_counts);
        rt_pdf(rt_pdf == 0) = nan;
        hr = entropyfun(rt_pdf);
        
        % compute contingency
        c = (hb - hr) / hb;
        
        % annotate entropy
        text(max(xlim)*1.05,offset,sprintf('C = %.2f',c),...
            'color',session_clrs(ss,:));
    end
    
    % adjust y-limits
    ylim(ylim + [-1,1] * .05 * range(ylim));

    % save figure
    if ~isempty(opt.savepath)
        file_name = sprintf('%s',[get(gcf,'name'),'.png']);
        file_path = fullfile(opt.savepath,file_name);
        print(gcf,file_path,'-dpng','-r300','-painters');
    end
end
end