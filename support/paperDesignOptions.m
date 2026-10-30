function opts = paperDesignOptions()
%PAPERDESIGNOPTIONS Shared dimensionless MB/DD design and experiment settings.
opts.TimeScale=1e-3; % normalized time tau=t/TimeScale
opts.StateScale=48; % same scale on both local states preserves scalar Xi blocks
opts.InputScale=48; % current and converter-command channels
opts.Rho=0.2;
opts.LocalGainLimit=20;
opts.GlobalGainLimit=20;
opts.Margin=1e-7;
opts.EdgePenalty=0.02;
opts.CommunicationBudget=10;
opts.PerformanceWeight=1;
opts.LocalRegularization=1e-3;
opts.AllowedEdges=[];
end
