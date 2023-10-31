function plotdopamineaverages_v2(data,meta,opt)
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

% iterate through mice
for mm = 1 : meta.mice.n
    mouse_flags = data.mouse == meta.mice.ids{mm};
    
    % parse session data
    n_sessions = max(data.session(mouse_flags));
    session_clrs = cool(n_sessions);
    
    % figure initialization
    figure(...
        'numbertitle','off',...
        'name',sprintf('%s_%s_da_averages_%s',...
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
    ylabel('DA (\DeltaF/F)');
    
    % mouse-specific DA range
    da_range = diff(...
        quantile(data.da.(opt.alignment)(mouse_flags,:),...
        [0,1]+[1,-1]*.05,'all'));
    
    % iterate through sessions
    for ss = 1 : n_sessions
        session_flags = data.session == ss;
        trial_flags = ...
            mouse_flags & ...
            session_flags;
        n_trials = sum(trial_flags);
        
        % vertical offset
        offset = -ss * da_range;
        
        % plot average DA
        da_mat = data.da.(opt.alignment)(trial_flags,:);
        da_time = linspace(...
            meta.epochs.(opt.alignment)(1),...
            meta.epochs.(opt.alignment)(2),size(da_mat,2));
        da_mu = nanmean(da_mat,1);
        da_std = nanstd(da_mat,0,1);
        da_sem = da_std ./ sqrt(n_trials);
        nan_flags = isnan(da_mu) | isnan(isnan(da_sem)) | da_sem == 0;
        errorpatch(da_time,da_mu+offset,da_sem,session_clrs(ss,:),...
            'edgecolor','none',...
            'facealpha',.25);
        plot(da_time(~nan_flags),da_mu(~nan_flags)+offset,...
            'color',session_clrs(ss,:),...
            'linewidth',1);
    
        % plot median reaction time
        rt_sign = (-1) ^ strcmpi(opt.alignment,'collection');
        rt_med = nanmedian(data.rt(trial_flags)) * rt_sign;
        idx = find(da_time >= min(rt_med, max(xlim)),1);
        if rt_med > max(xlim)
            marker = '>';
        else
            marker = 'v';
        end
        plot(da_time(idx),da_mu(idx)+offset,...
            'marker',marker,...
            'markersize',7.5,...
            'markeredgecolor','k',...
            'markerfacecolor',session_clrs(ss,:),...
            'linewidth',1.5);
    end
    
    % pseudo legend
    xleg = min(xlim)+.075*range(xlim);
    yleg = max(ylim)-.0*range(ylim);
    plot(xleg,yleg,...
        'color','none',...
        'marker','v',...
        'markersize',7.5,...
        'markeredgecolor','k',...
        'markerfacecolor','w',...
        'linewidth',1.5);
    text(xleg+.05*range(xlim),yleg,...
        '$\hat{rt}$',...
        'horizontalalignment','left',...
        'interpreter','latex',...
        'fontsize',12);
    
    % adjust y-limits
    ylim(ylim + [-1,1] * .05 * range(ylim));
    
    % plot alignment line
    plot([0,0],ylim,'--k');
    
    % save figure
    if ~isempty(opt.savepath)
        file_name = sprintf('%s',[get(gcf,'name'),'.png']);
        file_path = fullfile(opt.savepath,file_name);
        print(gcf,file_path,'-dpng','-r300','-painters');
    end
end
end