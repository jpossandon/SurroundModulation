function E354_checkTriggers(EEG,TRIAL)

%%
% numbers used as start trial triggers (>200 are target trials)
triggers = unique([TRIAL.trigger]);

first_trig_trial = 0;
firstTrigLatency = [];
ntriggers        = [];
misflickering    = [];
ntrg_pertrial    = 0;
inTrial          = 0;
for trg = 1:length(EEG.event)
    
    if ismember(str2num(EEG.event(trg).type(2:end)),triggers)
        first_trig_trial      = 1;
        if inTrial
            ntriggers         = [ntriggers ntrg_pertrial];
            if any(~ismember(unique(flickerLatencies),[33,34]))
                misflickering = [misflickering 1];
            else
                misflickering = [misflickering 0];
            end
            ntrg_pertrial    = 0;
            inTrial          = 0;
        end
        continue
    end
    if first_trig_trial ==1 & ismember(str2num(EEG.event(trg).type(2:end)),[30 80])
        firstTrigLatency = [firstTrigLatency,EEG.event(trg).latency-EEG.event(trg-1).latency];
        ntrg_pertrial    = ntrg_pertrial+1;
        inTrial          = 1;
        flickerLatencies = [];
        first_trig_trial = 0; 
    end
    if inTrial
        if ismember(str2num(EEG.event(trg).type(2:end)),[30 80])
            ntrg_pertrial    = ntrg_pertrial+1;
            flickerLatencies = [flickerLatencies,EEG.event(trg).latency-EEG.event(trg-1).latency];
        else
            ntriggers        = [ntriggers ntrg_pertrial ];
            if any(~ismember(unique(flickerLatencies),[33,34]))
                misflickering = [misflickering 1];
            else
                misflickering = [misflickering 0];
            end
            ntrg_pertrial    = 0;
            inTrial = 0;
        end
    end

end
