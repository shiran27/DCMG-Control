function plotLiteraturePIComparison(r,outdir)
%PLOTLITERATUREPICOMPARISON Compare PI and both co-designs on identical inputs.
indices=[2,3,r.LiteratureIndex];
labels={'MB co-design','DD local / known-B global','Passivity PI (Nahata et al.)'};
colors=[0.85,0.33,0.10;0.93,0.69,0.13;0.49,0.18,0.56];
fig=figure('Visible','off','Color','w','Position',[50,50,950,650]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for q=1:4
    ax=nexttile(tl,q); hold(ax,'on'); grid(ax,'on');
    for j=1:numel(indices)
        k=indices(j); style='-'; if j==3, style='--'; end
        if q==1
            c=r.Nominal{k}; y=vecnorm(c.StateError/48,2,2);
        elseif q==2
            c=r.Disturbance{k}; y=sqrt(mean(c.StateError(:,1:2:end).^2,2));
        elseif q==3
            c=r.Disturbance{k}; I=c.State(:,2:2:end)./r.PhysicalParameters.Irated;
            y=sqrt(mean((I-mean(I,2)).^2,2));
        else
            c=r.LoadStep{k}; y=sqrt(mean(c.StateError(:,1:2:end).^2,2));
        end
        if q==1 || q==4
            semilogy(ax,c.Time*1e3,max(y,1e-10),style,'Color',colors(j,:), ...
                'LineWidth',1.5,'DisplayName',labels{j}); set(ax,'YScale','log');
        else
            plot(ax,c.Time*1e3,y,style,'Color',colors(j,:), ...
                'LineWidth',1.5,'DisplayName',labels{j});
        end
    end
    xlim(ax,[0,120]); xlabel(ax,'Time (ms)');
    if q==1
        title(ax,'Nominal recovery'); ylabel(ax,'Normalized state-error norm'); ylim(ax,[1e-10,1]);
    elseif q==2
        title(ax,'Matched pulse and state jump'); ylabel(ax,'Voltage RMS error (V)');
    elseif q==3
        ylabel(ax,'Current-sharing RMS error (p.u.)');
    else
        title(ax,'Persistent load increment: +0.35 A per bus');
        ylabel(ax,'Voltage RMS error (V)'); ylim(ax,[1e-7,2]);
    end
    if q==2 || q==3, xline(ax,[25,35,65],':','HandleVisibility','off');
    elseif q==4, xline(ax,25,':','HandleVisibility','off'); end
end
lg=legend(ax,'Orientation','horizontal'); lg.Layout.Tile='south';
set(findall(fig,'Type','axes'),'FontSize',10,'Color','w','XColor','k','YColor','k');
set(findall(fig,'Type','text'),'Color','k');
name=fullfile(outdir,'literature_pi_comparison');
exportgraphics(fig,[name,'.pdf'],'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,[name,'.png'],'Resolution',220,'BackgroundColor','white');
savefig(fig,[name,'.fig']); close(fig);
end
