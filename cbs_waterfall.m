%% function waterfallplot

function [cfg] = cbs_waterfall(cfg,psd)

%make default plot
h = waterfall(psd);
ax = gca;

%change axis limits
if isfield(cfg,'xlim') && size(cfg.xlim,2) == 2;     ax.XLim = cfg.xlim;     end
if isfield(cfg,'ylim') && size(cfg.ylim,2) == 2;     ax.YLim = cfg.ylim;     end
if isfield(cfg,'zlim') && size(cfg.zlim,2) == 2;     ax.ZLim = cfg.zlim;     end

%change axis ticklabels
if isfield(cfg,'XTickLabel');     ax.XTickLabel = cfg.XTickLabel;     end
if isfield(cfg,'YTickLabel');     ax.YTickLabel = cfg.YTickLabel;     end
if isfield(cfg,'ZTickLabel');     ax.ZTickLabel = cfg.ZTickLabel;     end

%change axis labels
if isfield(cfg,'xlabel') && ischar(cfg.xlabel);     ax.XLabel.String = cfg.xlabel;  else    ax.XLabel.String = 'Frequency [Hz]';    end
if isfield(cfg,'ylabel') && ischar(cfg.ylabel);     ax.YLabel.String = cfg.ylabel;  else    ax.YLabel.String = 'Trials';            end
if isfield(cfg,'zlabel') && ischar(cfg.zlabel);     ax.ZLabel.String = cfg.zlabel;  else    ax.ZLabel.String = 'Power';             end

colormap('jet');
view(45,40);

h.LineWidth = 2;

end