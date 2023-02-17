% The script requires the Signal Processing Toolkit (Butter and filtfilt)
%The script requires the following functions: 
%abfload: https://github.com/fcollman/abfload
%triggerplot:https://github.com/VUntiet/IPSC-analysis/commit/30c325164773876d467707fb331d0ed6f097e3f0
%TiggerPoints:https://github.com/VUntiet/IPSC-analysis/commit/30c325164773876d467707fb331d0ed6f097e3f0
%TriggerPointsEnd: https://github.com/VUntiet/IPSC-analysis/commit/30c325164773876d467707fb331d0ed6f097e3f0

% zap all to make sure analysis is done anew.
clear; close all;

% ExperimentType definition: ONECOLOR NphR; TWOCOLOR=Swich
ONECOLOR=1;TWOCOLOR=2;NOSTIM=0;

% from CTN@KU
% DataDir = 'C:\Users\qzc716\Dropbox (CTN)\Matlab codes\Verena';
% from home
DataDir = 'H:\Work\Data\Rochester\IPSC for chloride paper\Verena';

DataFiles1 = {
    {  
    '11082022NphR31108\22n08003-ES-100uA-H-60mV.abf',ONECOLOR;
    }
    };

DataFiles10 = {
    {
    '12032022SwiChR3-1203\Slice3\22d03009-ES-100uA-H-40mV.abf',TWOCOLOR;
    }
    };



DataFiles=DataFiles1;


NumSlice = size(DataFiles,1);
for ii=1:NumSlice
    NumExperiment = size(DataFiles{ii},1);
    for jj=1:NumExperiment
        RawDataFileName=[DataDir,'\',DataFiles{ii}{jj,1}];
        [RawData,RawDataSampleInterval,RawDataInfo] = abfload(RawDataFileName);
        TAxis = [0:size(RawData,1)-1]'*RawDataSampleInterval*1E-6; %conversion to s (originally in us)
        SamplingFrequency = 1/RawDataSampleInterval * 1E6; %conversion to Hz;
        % FilteredEphysData = bandpass(RawData(:,1),[0.1,2000],SamplingFreqency);
        
        %preprocessing of raw data for ipsc detection
        [BB,AA] = butter(2,[0.5,300]/SamplingFrequency*2);
        EphysCap = -500; %pA
        CappedEphysData = RawData(:,1);
        CappedEphysData(CappedEphysData<EphysCap) = EphysCap;% EphysCap;
        FilteredEphysData = filtfilt(BB,AA,CappedEphysData);

        ExperimentType = DataFiles{ii}{jj,2};
        switch ExperimentType
            case ONECOLOR
                LightOn = TriggerPoints(RawData(:,3),1,100); % 1V threshold, 100 samples@10kHz = 10 ms minimum interval;
                LightOff = TriggerPointsEnd(RawData(:,3),1,100);

            case TWOCOLOR
                LightOn = TriggerPoints(RawData(:,2),1,100); % 1V threshold, 100 samples@10kHz = 10 ms minimum interval;
                LightOff = TriggerPoints(RawData(:,4),1,100);
            otherwise
                error('ExperimentType undefined')
        end

        APThreshold = 100; %pA if a peak is bigger than this, it is ignored. (i.e. action potential);
        MinPeakHeight = 10; %pA
        MinProminenceHeight = 5; %pA
        MinPeakDistanceMs = 50; %ms;
        MinPeakWidthMs = 2.5; %ms

        MinPeakDistanceSpl = MinPeakDistanceMs * 1E-3 * SamplingFrequency;
        MinPeakWidthSpl = MinPeakWidthMs * 1E-3 * SamplingFrequency;

        [Pks,Locs,Widths,Proms] = findpeaks(-FilteredEphysData,...
            'MinPeakHeight',MinPeakHeight,...
            'MinPeakDistance',MinPeakDistanceSpl,...
            'MinPeakProminence',MinProminenceHeight,...
            'MinPeakWidth',MinPeakWidthSpl);

        NonAPLocs=Locs(Pks<APThreshold);
        NonAPWidths=Widths(Pks<APThreshold);
        NonAPProms=Proms(Pks<APThreshold);

        Off1Locs = NonAPLocs(NonAPLocs < LightOn);
        Off1Widths = NonAPWidths(NonAPLocs < LightOn)./SamplingFrequency*1E3; %ms;
        Off1Proms = NonAPProms(NonAPLocs < LightOn);
        OnLocs = NonAPLocs((NonAPLocs >= LightOn) & (NonAPLocs < LightOff));
        OnWidths = NonAPWidths((NonAPLocs >= LightOn) & (NonAPLocs < LightOff))./SamplingFrequency*1E3; %ms;
        OnProms = NonAPProms((NonAPLocs >= LightOn) & (NonAPLocs < LightOff));
        Off2Locs = NonAPLocs(NonAPLocs >= LightOff);
        Off2Widths = NonAPWidths(NonAPLocs >= LightOff)./SamplingFrequency*1E3; %ms;
        Off2Proms = NonAPProms(NonAPLocs >= LightOff);

        tmpfig = figure;
        ax1= subplot(2,2,1);
        plot(TAxis,FilteredEphysData); hold on;
        plot(TAxis(Off1Locs),FilteredEphysData(Off1Locs),'ro');
        plot(TAxis(OnLocs),FilteredEphysData(OnLocs),'go');
        plot(TAxis(Off2Locs),FilteredEphysData(Off2Locs),'mo');
        xlabel ('time (s)'); ylabel ('pA');

        ax2= subplot(2,2,3);
        plot(TAxis,detrend(CappedEphysData,'constant')); hold on;
        plot(TAxis,CappedEphysData,'r'); hold on;
        xlabel ('time (s)'); ylabel ('pA');

        subplot(2,2,2)
        %Data4TriggerPlot = RawData(:,1);
        Data4TriggerPlot = FilteredEphysData(:,1);
        [Off1TracesMean,Off1TracesSTD,t,Off1TracesAll] = triggerplot(Data4TriggerPlot,Off1Locs,500,1500);
        [OnTracesMean,OnTracesSTD,t,OnTracesAll] = triggerplot(Data4TriggerPlot,OnLocs,500,1500);
        [Off2TracesMean,Off2TracesSTD,t,Off2TracesAll] = triggerplot(Data4TriggerPlot,Off2Locs,500,1500);
        TriggeredTAxisMs = t./SamplingFrequency*1E3; %ms
        plot(TriggeredTAxisMs,Off1TracesMean,'r'); hold on;
        plot(TriggeredTAxisMs,OnTracesMean,'g'); hold on;
        plot(TriggeredTAxisMs,Off2TracesMean,'m'); hold on;
        xlabel ('time (ms)'); ylabel ('pA');

        linkaxes([ax1,ax2]);

        % saveas(tmpfig,sprintf('slice-%d-trial-%d.png',ii,jj));
          % modified 20220917 to plot cumulative distributions%
        BeforeLength = 500; % 500 points@10 kHz = 50 ms.
        AfterLength = 1500; % 1500 points@10 kHz = 150 ms.
        [Off1Traces,s,t,Off1TraceWave] = triggerplot(Data4TriggerPlot,Off1Locs,BeforeLength,AfterLength);
        [OnTraces,~,~,OnTraceWave] = triggerplot(Data4TriggerPlot,OnLocs,BeforeLength,AfterLength);
        [Off2Traces,~,~,Off2TraceWave] = triggerplot(Data4TriggerPlot,Off2Locs,BeforeLength,AfterLength);
        TriggeredTAxisMs = t ./SamplingFrequency*1E3; %ms
        
        plot(TriggeredTAxisMs,Off1Traces,'r'); hold on;
        plot(TriggeredTAxisMs,OnTraces,'g'); hold on;
        plot(TriggeredTAxisMs,Off2Traces,'m'); hold on;
        xlabel ('time (ms)'); ylabel ('pA');

        linkaxes([ax1,ax2]);

        %SAVES DATA
        [filepath,name,ext] = fileparts(RawDataFileName); % Get name of file to save tables

        % added 20220917 draw cumulative distributions (needs stats toolbox)
        subplot(2,2,4)
        % peak is detected at BeforeLength+1 in the first dimension of *Wave.
        [Off1CDF,Off1CDFX] = ecdf(abs(Off1TraceWave(BeforeLength+1,:)));
        [OnCDF,OnCDFX] = ecdf(abs(OnTraceWave(BeforeLength+1,:)));
        [Off2CDF,Off2CDFX] = ecdf(abs(Off2TraceWave(BeforeLength+1,:)));
         tab_Off1CDF_params = array2table([Off1CDF,Off1CDFX],'VariableNames',{'Off1CDFCDF (s)', 'Off1CDFX (s)'});
        writetable(tab_Off1CDF_params,[filepath filesep name '_Off1CDF_params_All']);
          tab_Off2CDF_params = array2table([Off2CDF,Off2CDFX],'VariableNames',{'Off2CDFCDF (s)', 'Off2CDFX (s)'});
        writetable(tab_Off2CDF_params,[filepath filesep name '_Off2CDF_params_All']);
          tab_OnCDF_params = array2table([OnCDF,OnCDFX],'VariableNames',{'OnCDFCDF (s)', 'OnCDFX (s)'});
        writetable(tab_OnCDF_params,[filepath filesep name '_OnCDF_params_All']);
        plot(Off1CDFX,Off1CDF,'r'); hold on;
        plot(OnCDFX,OnCDF,'g'); hold on;
        plot(Off2CDFX,Off2CDF,'m'); hold on;
        % perform rank sum test between Pre (red) and During (green).
        p = ranksum(abs(Off1TraceWave(BeforeLength+1,:)),abs(OnTraceWave(BeforeLength+1,:)));
        xlim([15,100]); xlabel('IPSC (pA)'); ylabel('Cumulative probab.');
        title(sprintf('pre vs. during p=%5.3f',p));
        % saveas(tmpfig,sprintf('slice-%d-trial-%d.png',ii,jj));
        
        % OFF 1 
        % Mean Trace
        tab_Off1_mean = array2table([TriggeredTAxisMs' Off1TracesMean Off1TracesSTD],'VariableNAmes',{'Time (ms)' 'IPSC Mean' 'IPSC STD'});
        writetable(tab_Off1_mean,[filepath filesep name '_Off1_Traces_Mean']);
        % All Traces
        tab_Off1_all = array2table([TriggeredTAxisMs' Off1TracesAll]);
        tab_Off1_all = renamevars(tab_Off1_all,'Var1','Time (ms)');
        writetable(tab_Off1_all,[filepath filesep name '_Off1_Traces_All']);
        % Parameters: Amplitude, Frequency, Widths, Proms, Rise&Decay Time
        Off1_Traces_Amplitude_All = FilteredEphysData(Off1Locs);
        Off1_Traces_Frequency_All = 1./[diff(Off1Locs)./SamplingFrequency]; %Hz
        Off1_Traces_Frequency_All = [Off1_Traces_Frequency_All; Off1_Traces_Frequency_All(end)]; %Hz
        Off1_Traces_Widths_All = Off1Widths;
        Off1_Traces_Prominence_All = Off1Proms;
        [Off1_Traces_TRise_All, Off1_Traces_TDecay_All] = rise_decay_time(FilteredEphysData, Off1Locs, SamplingFrequency);
        tab_Off1_params = array2table([Off1Locs/SamplingFrequency,Off1_Traces_Amplitude_All, Off1_Traces_Frequency_All, Off1_Traces_Widths_All, Off1_Traces_Prominence_All, Off1_Traces_TRise_All, Off1_Traces_TDecay_All],'VariableNames',{'Peak Time (s)', 'Amplitude (pA)', 'Frequency (Hz)', 'Width (ms)', 'Prominence (pA)', 'Rise Time (ms)', 'Decay Time (ms)'});
        writetable(tab_Off1_params,[filepath filesep name '_Off1_Traces_Parameters_All']);
                
        % ON
        % Mean Trace
        tab_On_mean = array2table([TriggeredTAxisMs' OnTracesMean OnTracesSTD],'VariableNAmes',{'Time (ms)' 'IPSC Mean' 'IPSC STD'});
        writetable(tab_On_mean,[filepath filesep name '_On_Traces_Mean']);
        % All Traces
        tab_On_all = array2table([TriggeredTAxisMs' OnTracesAll]);
        tab_On_all = renamevars(tab_On_all,'Var1','Time (ms)');
        writetable(tab_On_all,[filepath filesep name '_On_Traces_All']);
        % Parameters: Amplitude, Frequency, Widths, Proms, Rise&Decay Time
        On_Traces_Amplitude_All = FilteredEphysData(OnLocs);
        On_Traces_Frequency_All = 1./[diff(OnLocs)./SamplingFrequency]; %Hz
        On_Traces_Frequency_All = [On_Traces_Frequency_All; On_Traces_Frequency_All(end)]; %Hz
        On_Traces_Widths_All = OnWidths;
        On_Traces_Prominence_All = OnProms;
        [On_Traces_TRise_All, On_Traces_TDecay_All] = rise_decay_time(FilteredEphysData, OnLocs, SamplingFrequency);
        tab_On_params = array2table([OnLocs/SamplingFrequency,On_Traces_Amplitude_All, On_Traces_Frequency_All, On_Traces_Widths_All, On_Traces_Prominence_All, On_Traces_TRise_All, On_Traces_TDecay_All],'VariableNames',{'Peak Time (s)', 'Amplitude (pA)', 'Frequency (Hz)', 'Width (ms)', 'Prominence (pA)', 'Rise Time (ms)', 'Decay Time (ms)'});
        writetable(tab_On_params,[filepath filesep name '_On_Traces_Parameters_All']); 
        
        % OFF 2
        % Mean Trace
        tab_Off2_mean = array2table([TriggeredTAxisMs' Off2TracesMean Off2TracesSTD],'VariableNAmes',{'Time (ms)' 'IPSC Mean' 'IPSC STD'});
        writetable(tab_Off2_mean,[filepath filesep name '_Off2_Traces_Mean']);
        % All Traces
        tab_Off2_all = array2table([TriggeredTAxisMs' Off2TracesAll]);
        tab_Off2_all = renamevars(tab_Off2_all,'Var1','Time (ms)');
        writetable(tab_Off2_all,[filepath filesep name '_Off2_Traces_All']);
        % Parameters: Amplitude, Frequency, Widths, Proms, Rise&Decay Time
        Off2_Traces_Amplitude_All = FilteredEphysData(Off2Locs);
        Off2_Traces_Frequency_All = 1./[diff(Off2Locs)./SamplingFrequency]; %Hz
        Off2_Traces_Frequency_All = [Off2_Traces_Frequency_All; Off2_Traces_Frequency_All(end)]; %Hz
        Off2_Traces_Widths_All = Off2Widths;
        Off2_Traces_Prominence_All = Off2Proms;
        [Off2_Traces_TRise_All, Off2_Traces_TDecay_All] = rise_decay_time(FilteredEphysData, Off2Locs, SamplingFrequency);
        tab_Off2_params = array2table([Off2Locs/SamplingFrequency,Off2_Traces_Amplitude_All, Off2_Traces_Frequency_All, Off2_Traces_Widths_All, Off2_Traces_Prominence_All, Off2_Traces_TRise_All, Off2_Traces_TDecay_All],'VariableNames',{'Peak Time (s)', 'Amplitude (pA)', 'Frequency (Hz)', 'Width (ms)', 'Prominence (pA)', 'Rise Time (ms)', 'Decay Time (ms)'});
        writetable(tab_Off2_params,[filepath filesep name '_Off2_Traces_Parameters_All']);
              
        % Histogram of Amplitudes
        figure
        ax1=subplot(1,3,1);
        histogram(FilteredEphysData(Off1Locs),'FaceColor','r')
        grid on
        xlabel('Peak Amplitude (pA)')
        ylabel('Occurrence')
        title('Off 1 Interval')
        ax2=subplot(1,3,2);
        histogram(FilteredEphysData(OnLocs),'FaceColor','g')
        grid on
        xlabel('Peak Amplitude (pA)')
        ylabel('Occurrence')
        title('ON Interval')
        ax3=subplot(1,3,3);
        histogram(FilteredEphysData(Off2Locs),'FaceColor','m')
        grid on
        xlabel('Peak Amplitude (pA)')
        ylabel('Occurrence')
        title('Off 2 Interval')
        linkaxes([ax1 ax2 ax3],'xy')    
        sgtitle('Histogram of Peak Amplitudes (pA)')

        
        %% Electrical Stimulation Analysis
        if (size(RawData,2)>4)&&(~isempty(find(RawData(:,5)>4))) % Check if there was electrical stimulation            
            [BB,AA] = butter(2,[50]/(SamplingFrequency/2));
            Data4TriggerPlot_ES = filtfilt(BB,AA,RawData(:,1));
            
            % ES 1 -> 5 - 20 s
            ind_es1 = find(TAxis>5 & TAxis<20);
            ES1On = TriggerPoints(RawData(ind_es1,5),1,100); % 4V threshold, ES@50Hz (20ms interval @ 10KHz = 200samples);
            ES1On = ind_es1(ES1On);
            ES1Off = TriggerPointsEnd(RawData(ind_es1,5),1,100);
            ES1Off = ind_es1(ES1Off);
            
            % ES 2 -> 60 - 80 s
            ind_es2 = find(TAxis>60 & TAxis<80);
            ES2On = TriggerPoints(RawData(ind_es2,5),1,100); % 4V threshold, ES@50Hz (20ms interval @ 10KHz = 200samples);
            ES2On = ind_es2(ES2On);
            ES2Off = TriggerPointsEnd(RawData(ind_es2,5),1,100);
            ES2Off = ind_es2(ES2Off);
                        
            % ES 3 -> 100 - 120 s
            ind_es3 = find(TAxis>100 & TAxis<120);
            ES3On = TriggerPoints(RawData(ind_es3,5),1,100); % 4V threshold, ES@50Hz (20ms interval @ 10KHz = 200samples);
            ES3On = ind_es3(ES3On);
            ES3Off = TriggerPointsEnd(RawData(ind_es3,5),1,100);
            ES3Off = ind_es3(ES3Off);
            
            % ES 4 -> 165 - 185 s
            ind_es4 = find(TAxis>165 & TAxis<185);
            ES4On = TriggerPoints(RawData(ind_es4,5),1,100); % 4V threshold, ES@50Hz (20ms interval @ 10KHz = 200samples);
            ES4On = ind_es4(ES4On);
            ES4Off = TriggerPointsEnd(RawData(ind_es4,5),1,100);
            ES4Off = ind_es4(ES4Off);
            
            % Get the amplitude of the tracer at each stimulation of the burst 
            ESNumber = 1:100; % There are 100 ES per burst.
            ES1Locs = round(mean([ES1On ES1Off],2)); 
            ES2Locs = round(mean([ES2On ES2Off],2));
            ES3Locs = round(mean([ES3On ES3Off],2));
            ES4Locs = round(mean([ES4On ES4Off],2));
            ES1Amplitudes = Data4TriggerPlot_ES(ES1Locs);
            ES2Amplitudes = Data4TriggerPlot_ES(ES2Locs);
            ES3Amplitudes = Data4TriggerPlot_ES(ES3Locs);
            ES4Amplitudes = Data4TriggerPlot_ES(ES4Locs);
            
            % Plot Tracer during the 4 different stimulation bursts.
            figure
            subplot(3,1,1)
            [ES1Trace,~,t,~] = triggerplot(Data4TriggerPlot_ES,ES1Locs(1),10000,30000);
            [ES2Trace,~,t,~] = triggerplot(Data4TriggerPlot_ES,ES2Locs(1),10000,30000);
            [ES3Trace,~,t,~] = triggerplot(Data4TriggerPlot_ES,ES3Locs(1),10000,30000);
            [ES4Trace,~,t,~] = triggerplot(Data4TriggerPlot_ES,ES4Locs(1),10000,30000);
            ESTAxisMs = t./SamplingFrequency*1E3; %ms
            c_lines = colormap(lines);
            plot(ESTAxisMs,ES1Trace,'Color',c_lines(1,:)); hold on;
            plot(ESTAxisMs,ES2Trace,'Color',c_lines(2,:)); hold on;
            plot(ESTAxisMs,ES3Trace,'Color',c_lines(3,:)); hold on;
            plot(ESTAxisMs,ES4Trace,'Color',c_lines(4,:)); hold on;
            xlabel ('Iime (ms)'); ylabel ('pA');
            xline(0)
            xline((ES1Locs(2)-ES1Locs(1))/SamplingFrequency*1E3,'--k') % Second Stimulus (200us)
            xline((ES1Locs(90)-ES1Locs(1))/SamplingFrequency*1E3,'--k') % 90th Stimulus (200us)
            xline(2000)
            hold off
            legend('1st burst', '2nd burst', '3rd burst', '4th burst')
            
            subplot(3,1,2)
            [ES1Stim,~,t,~] = triggerplot(RawData(:,5),ES1On(1),10000,30000);
            plot(ESTAxisMs,ES1Stim,'k')
            
            subplot(3,1,3)
            bar([ES1Amplitudes(2) ES2Amplitudes(2) ES3Amplitudes(2) ES4Amplitudes(2);ES1Amplitudes(90) ES2Amplitudes(90) ES3Amplitudes(90) ES4Amplitudes(90)]);
            row1 = {'Early phase' 'Late phase'};
            row2 = {'(2nd stimulus)' '(90th stimulus)'};
            labelArray = [row1; row2];
            tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
            ax = gca(); 
            ax.XTick = 1:2; 
            ax.XLim = [0,3];
            ax.XTickLabel = tickLabels;
            ylabel('pA')
            legend('1st burst', '2nd burst', '3rd burst', '4th burst')

            
            % Save the traces
            tab_ES_trace_mean = array2table([ESTAxisMs' ES1Trace ES2Trace ES3Trace ES4Trace],'VariableNAmes',{'Time (ms)' 'ES 1 Trace' 'ES 2 Trace' 'ES 3 Trace' 'ES 4 Trace'});
            writetable(tab_ES_trace_mean,[filepath filesep name '_ES_Traces_Mean']);
            % Save the amplitudes
            tab_ES_amplitude_mean = array2table([ESNumber' ES1Amplitudes ES2Amplitudes ES3Amplitudes ES4Amplitudes],'VariableNAmes',{'ES Number' 'ES 1 Amplitudes' 'ES 2 Amplitudes' 'ES 3 Amplitudes' 'ES 4 Amplitudes'});
            writetable(tab_ES_amplitude_mean,[filepath filesep name '_ES_Traces_Amplitude']);
            
        end

        
    end
end
    


tab_ES_trace_mean_ds = downsample(tab_ES_trace_mean,100);
FilteredEphysData_ds = downsample(FilteredEphysData,10);