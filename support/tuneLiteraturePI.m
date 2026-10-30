function [controller,tuning] = tuneLiteraturePI(net,initialState,peakCommandBudget)
%TUNELITERATUREPI Tune on nominal recovery only, before disturbance evaluation.
% The peak-command budget equals the larger accepted MB/DD nominal peak.
% Select fastest 2% settling, then smallest state energy, then command peak.
kp=[0.25,0.5,1,2,4]; kd=[0.25,0.5,1,2,4]; fractions=[0.1,0.25,0.5,0.75,0.9];
rows=zeros(numel(kp)*numel(kd)*numel(fractions),9); index=0;
best=[Inf,Inf,Inf]; controller=[];
for a=kp
    for b=kd
        for f=fractions
            index=index+1;
            c=designLiteraturePI(net,VoltageGain=a,CurrentGain=b,IntegralFraction=f);
            accepted=false; settle=NaN; energy=NaN; peak=NaN; final=NaN;
            if c.problem==0
                sim=simulatePaperCase(net,c,InitialState=initialState);
                m=sim.Metrics; settle=m.SettlingTime2Percent; energy=m.StateEnergy;
                peak=m.PeakCommandDeviation; final=m.FinalRatio;
                accepted=sim.Accepted && isfinite(settle) && final<1e-5 ...
                    && peak<=peakCommandBudget*(1+1e-10);
                score=[settle,energy,peak];
                first=find(score~=best,1);
                if accepted && ~isempty(first) && score(first)<best(first)
                    controller=c; best=score;
                end
            end
            rows(index,:)=[a,b,f,c.GainNorm,settle,energy,peak,final,accepted];
        end
    end
end
assert(~isempty(controller),'DCMG:PINoTuning','No PI candidate passed the declared budget.');
tuning=struct('Candidates',array2table(rows,'VariableNames',{'VoltageGain','CurrentGain', ...
    'IntegralFraction','NormalizedGainNorm','SettlingTime','StateEnergy', ...
    'PeakCommandDeviation','FinalErrorRatio','Accepted'}), ...
    'PeakCommandBudget',peakCommandBudget,'Selected',controller.Options, ...
    'SelectionRule','Nominal settling time, then state energy, then command peak; no disturbance tuning');
fprintf('Literature PI: kp=%.3g, kd=%.3g, fraction=%.3g, settle=%.3g ms (%d/%d accepted).\n', ...
    controller.Options.VoltageGain,controller.Options.CurrentGain, ...
    controller.Options.IntegralFraction,1e3*best(1),sum(rows(:,end)),index);
end
