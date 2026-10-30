function plotPaperResults(report,outdir)
%PLOTPAPERRESULTS Recreate all publication figures from the saved result struct.
% Can be run alone after load('Results/Regenerated/paper_results.mat').
set(groot,'defaultFigureColor','w'); colors=lines(6);
labels=cellstr(report.Methods);
numberOfMethods=numel(labels);
fig=figure('Visible','off','Position',[50,50,max(1200,360*numberOfMethods),670]);
tl=tiledlayout(fig,2,numberOfMethods,'TileSpacing','compact','Padding','compact');
for k=1:numberOfMethods
    c=report.ArchivedInitialState{k};
    ax=nexttile(tl,k); hold(ax,'on');
    for i=1:6
        plot(ax,1e3*c.Time,c.State(:,2*i-1)/c.SteadyState.x_ss(2*i-1), ...
            'Color',colors(i,:),'LineWidth',1.1);
    end
    yline(ax,1,'k:'); grid(ax,'on'); title(ax,labels{k},'Interpreter','none');
    xlabel(ax,'Time (ms)'); ylabel(ax,'Voltage / equilibrium voltage'); ylim(ax,[-1.5,2.5]);
    ax=nexttile(tl,numberOfMethods+k); hold(ax,'on');
    for i=1:6
        plot(ax,1e3*c.Time,c.State(:,2*i)/report.PhysicalParameters.Irated(i), ...
            'Color',colors(i,:),'LineWidth',1.1);
    end
    yline(ax,c.SteadyState.epsI,'k:'); grid(ax,'on');
    xlabel(ax,'Time (ms)'); ylabel(ax,'Converter current / rated current'); ylim(ax,[-4,4]);
end
lg=legend(ax,compose('DG %d',1:6),'Orientation','horizontal'); lg.Layout.Tile='south';
savePlot(fig,outdir,'paper_trajectories');

fig=figure('Visible','off','Position',[50,50,max(1200,360*numberOfMethods),450]);
tl=tiledlayout(fig,1,numberOfMethods,'TileSpacing','compact','Padding','compact');
for k=1:numberOfMethods
    [net,~]=createPaperDCMicrogrid(LegacyEquilibrium=false);
    net.K=report.Controllers{k}.CommandGain;
    net.setStateVector(net.x_s);
    % Draw the existing physical topology without crowded electrical labels.
    ax=nexttile(tl,k); hold(ax,'on'); p=report.PhysicalParameters.Positions;
    for j=1:net.M
        line=net.Lines(j);
        plot(ax,p([line.i,line.j],1),p([line.i,line.j],2),'-','Color',[.55,.55,.55],'LineWidth',2);
    end
    K=report.Controllers{k}.CommandGain; adj=false(6);
    for i=1:6
        for j=1:6, adj(i,j)=i~=j && norm(K(i,2*j-1:2*j))>1e-9; end
    end
    G=digraph(adj');
    plot(ax,G,'XData',p(:,1),'YData',p(:,2),'NodeColor',[0.05,.32,.48], ...
        'EdgeColor',[.9,.35,.08],'LineStyle','--','ArrowSize',7,'MarkerSize',9, ...
        'NodeLabel',compose('%d',1:6),'LineWidth',0.8);
    axis(ax,'equal'); axis(ax,'off');
    xlim(ax,[min(p(:,1))-.6,max(p(:,1))+.6]);
    ylim(ax,[min(p(:,2))-.6,max(p(:,2))+.6]);
    title(ax,sprintf('%s\n%d directed communication links',labels{k},nnz(adj)), ...
        'Interpreter','none');
end
title(tl,'Physical lines: gray solid     Communication: orange dashed','FontWeight','normal');
savePlot(fig,outdir,'paper_topologies');

fig=figure('Visible','off','Position',[50,50,1000,700]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for q=1:4
    ax=nexttile(tl,q); hold(ax,'on'); grid(ax,'on');
    for k=1:numberOfMethods
        if q==1
            c=report.Nominal{k}; y=vecnorm(c.StateError/48,2,2);
            semilogy(ax,1e3*c.Time,max(y,1e-12),'LineWidth',1.4,'DisplayName',labels{k});
        elseif q==2
            c=report.Disturbance{k}; y=sqrt(mean(c.StateError(:,1:2:end).^2,2));
            plot(ax,1e3*c.Time,y,'LineWidth',1.4,'DisplayName',labels{k});
        elseif q==3
            c=report.Disturbance{k}; I=c.State(:,2:2:end)./report.PhysicalParameters.Irated;
            y=sqrt(mean((I-mean(I,2)).^2,2));
            plot(ax,1e3*c.Time,y,'LineWidth',1.4,'DisplayName',labels{k});
        else
            c=report.Disturbance{k}; y=sqrt(mean((c.Command-c.SteadyState.u_s').^2,2));
            plot(ax,1e3*c.Time,y,'LineWidth',1.4,'DisplayName',labels{k});
        end
    end
    xlabel(ax,'Time (ms)'); xlim(ax,[0,120]);
    if q==1, ylabel(ax,'Normalized state-error norm'); set(ax,'YScale','log'); title(ax,'Nominal recovery');
    elseif q==2, ylabel(ax,'Voltage RMS error (V)'); title(ax,'Matched disturbance experiment');
    elseif q==3, ylabel(ax,'Current-sharing RMS error (p.u.)');
    else, ylabel(ax,'Command RMS deviation (V)'); end
    if q>1, xline(ax,[25,35,65],':','Color',[.6,.6,.6],'HandleVisibility','off'); end
end
lg=legend(ax,'Location','southoutside','Orientation','horizontal'); lg.Layout.Tile='south';
savePlot(fig,outdir,'paper_comparison');

if isfield(report,'Variation')
    fig=figure('Visible','off','Position',[50,50,900,360]);
    tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    for k=1:2
        ax=nexttile(tl,k); hold(ax,'on'); grid(ax,'on');
        f=report.Variation.Frozen{k+1}; r=report.Variation.Redesigned{k};
        plot(ax,1e3*f.Time,sqrt(mean(f.StateError(:,1:2:end).^2,2)),'LineWidth',1.3);
        plot(ax,1e3*r.Time,sqrt(mean(r.StateError(:,1:2:end).^2,2)),'LineWidth',1.3);
        xlabel(ax,'Time (ms)'); xlim(ax,[0,120]); ylabel(ax,'Voltage RMS error (V)');
        title(ax,labels{k+1},'Interpreter','none'); legend(ax,{'Frozen gains','Offline redesigned'},'Location','best');
    end
    title(tl,'Physical parameter variation within +/-10% (empirical stress test)');
    savePlot(fig,outdir,'paper_variation');
end
end

function savePlot(fig,outdir,name)
set(findall(fig,'Type','axes'),'FontSize',10,'Color','w','XColor','k','YColor','k');
set(findall(fig,'Type','text'),'Color','k');
exportgraphics(fig,fullfile(outdir,[name,'.png']),'Resolution',220,'BackgroundColor','white');
exportgraphics(fig,fullfile(outdir,[name,'.pdf']),'ContentType','vector','BackgroundColor','white');
savefig(fig,fullfile(outdir,[name,'.fig'])); close(fig);
end
