function [peak_height,std_value, peak_positions, activation, deactivation] = calculatePeakStatistics (data, mask, Options)

%Function used by script_peak_improved to find the peaks of the maxima
%signal

% V 1.0 Creation of the document by David Lopez 30.04.2015
% V 1.1 The New Dicom data has been added by David Lopez 30.07.2015
% V 1.2 The activation values have been added by David Lopez 09.08.2015
% V 1.3 If the array of pos is given for the 2D case then take those
% positions all the time By David Lopez 28.08.2015

%Function Inputs
%-data --> data.images for 3D or data.signal for 2D analysis
%-mask --> mask with the pixels for the 3D data
%-Options:
%+start_stim --> first point of the stimulous --> Future adaptation
%+dur_stim --> duration of the stimulation --> Future adaptation
%+data_type -> type of data to know the format

%Functions Outputs
%- peak_height --> returns the 2D map with the average peak height for each
%voxel or the value for the individual series
%- std_value --> returns the 2D map with the std peak height for each
%voxel or the value for the individual series
%- peak_height_act --> values during the activation period --> Future adaptation
%- peak_height_deact --> values during the rest period --> Future adaptation
%- peak_positions -> returns the array with the positions of the peak
%maxima


%Check the initial parameters
if nargin < 1
    error('Signal or data arguments need to be provided');
end

if (nargin < 2) && isfield(data,'images')
    error('A mask is need to carried out the slice analysis');
end

if nargin < 3
    activation = [];
    deactivation = [];
    activationFlag = 0;
elseif (nargin == 3)
    if isfield(Options,'start_stim') &&  ~(isfield(Options,'dur_stim'))
        error('The beginning of the stimulous has been set but not the duration');
    elseif ~isfield(Options,'start_stim') &&  (isfield(Options,'dur_stim'))
        error('The beginning of the stimulous has not been set but the duration has');
    elseif isfield(Options,'start_stim') &&  (isfield(Options,'dur_stim'))
        activationFlag = 1;
    else
        activation = [];
        deactivation = [];
        activationFlag = 0;
    end
end

if isfield(data,'images')
    if isfield(Options,'positions')
        case_type = '2DPeaks';
    else
        case_type = '2D';
    end
elseif isfield(data,'signal')
    case_type = 'Signal';
else
    error('The signal or the data has to be provided');
end

%Common variables
if strcmp(Options.DataType, 'Philips')
    dt = str2num(data.parameters.Repetitiontimems)/1000;% Time spacing
elseif strcmp(Options.DataType, 'Dicom')
    dt = data.parameters.RepetitionTime/1000;% Time spacing
else
    if strcmp(case_type,'2D') || strcmp(case_type,'2DPeaks')
        dt = data.parameters.AcquisitionDuration/size(data.images,3);% Time spacing   
    else
        dt = data.parameters.AcquisitionDuration/size(data.signal,2);% Time spacing   
    end
end
Fs = 1/dt;% Frequency spacing
    
%Process
switch case_type 
    %Slice case
    case ('2D')
        cardiac_peak_map = zeros(size(mask));
        N = size(data.images,3);% Number of files
        dF = Fs/N;% Frequency spacing
        f = -Fs/2:dF:Fs/2-dF;% Frequency scale  
        BPF = ((0.6 < abs(f)) & (abs(f) < 1.5));% To detect the cardiac maxima
        peak_height = zeros(size(mask));%Initialize mean array
        std_value = zeros(size(mask));%Initialize std array
        
        for i=1:size(data.images,1)
            for j=1:size(data.images,2)
                if mask(i,j)~=0
                    signal_cardiac(1,:) =  double(data.images(i,j,:));
                    signal_cardiacaux(1,:)=detrend(signal_cardiac,'constant');      
                    spektrum = ifftshift (abs(fft(signal_cardiacaux)));
                    spektrum = BPF.*spektrum;
                    [x,y] = find(spektrum==max(spektrum));
                    cardiac_peak_map(i,j) = spektrum(1,y(1,1));
                    %Divide the signal into blocks and calculate the maxima of
                    %the block
                    cardiac_freq = dF*(size(f,2)/2-y(1,1)+1);
                    OptionsPeak.block_size = round((1/cardiac_freq)/dt);
                    blocks = round(size(signal_cardiac,2)/OptionsPeak.block_size);
                    OptionsPeak.ini = 0;
                    OptionsPeak.estimatedInterval = round(300/(dt*1000))+1;%Variations up to 200ms in time
                    mean_signal = mean(signal_cardiacaux(signal_cardiacaux>0));
                    signal_interval = signal_cardiacaux((1):(1+OptionsPeak.block_size));
                    %Location of the first peak
                    [pos(1),~,ini] = peakdet(signal_interval, 1.5*mean_signal,OptionsPeak);                 
                    if ini <= OptionsPeak.estimatedInterval
                        OptionsPeak.ini = ini+OptionsPeak.block_size;
                        pos(1) = 0;
                    else
                        OptionsPeak.ini = ini;
                    end
                    clear signal_interval;
                    %Once the first peak is located search for the rest of the peaks
                    OptionsPeak.estimatedInterval = round(300/(dt*1000));%Variations up to 200ms in time
                    for b = 1:blocks-5% The last 5 blocks are discarded to avoid problems with the array length
                        if ((OptionsPeak.ini+OptionsPeak.estimatedInterval)<size(signal_cardiac,2))
                            signal_interval = signal_cardiac((OptionsPeak.ini-OptionsPeak.estimatedInterval):(OptionsPeak.ini+OptionsPeak.estimatedInterval));
                            [pos(b+1),~,OptionsPeak.ini] = peakdet(signal_interval,1.5*mean_signal,OptionsPeak);
                        end
                    end
                    value_total = 0;
                    pos(pos==0) = [];%Removed positions where no peak was found
                    %Calculate the total mean peak
                    for k = 1:size(pos,2)        
                        value_total = double(value_total) + signal_cardiac(1,pos(k));
                    end

                    mean_signal2 = double(mean(signal_cardiac(signal_cardiac>0)));%Mean value of the whole signal
                    mean_peak = double(value_total)/size(pos,2);%Mean value for all the peaks

                    approximateMean = (double(mean_signal2)*(size(signal_cardiac,2)/(N-size(pos,2))))- (double(value_total)/(N-size(pos,2)));%Aprox. mean without peaks
                    peak_height(i,j) = double(mean_peak) - double(approximateMean);
                    std_value(i,j) = std(signal_cardiac(pos));%Std of the peaks
                    peak_positions(i,j).pos = pos;
                    %Classify the peaks among activation and rest
                    if activationFlag
                        numberOfBlocks = round((size(data.images,3)-Options.start_stim)/Options.dur_stim)/2;%number of activation/rest blocks
                        array_aux = dt:dt:(Options.dur_stim*dt);
                        posAct = zeros(size(pos));
                        posRest = zeros(size(pos));
                        for nb = 1:numberOfBlocks
                            %Activation
                            activationRange = ((Options.start_stim + Options.dur_stim*((nb-1)*2))+1) : (Options.start_stim + Options.dur_stim*(((nb-1)*2)+1));
                            %Fitting of the range
                            signalFit(1,:) = data.images(i,j,activationRange);
                            [pAct, rAct] = polyfit(array_aux,double(signalFit),1); 
                            %Rest
                            restRange = ((Options.start_stim + Options.dur_stim*(((nb-1)*2)+1))+1) : (Options.start_stim + Options.dur_stim*(((nb)*2)));
                            signalFit(1,:) = data.images(i,j,restRange);
                            [pRest, rRest] = polyfit(array_aux,double(signalFit),1); 
                            %Check the pos array for those values in the
                            %activation or rest range
                            for actresPos = 1:size(pos,2)
                                %if act
                                if ((pos(actresPos) >= (Options.start_stim + Options.dur_stim*((nb-1)*2))+1) && (pos(actresPos) <= (Options.start_stim + Options.dur_stim*(((nb-1)*2)+1))))
                                    posAct(actresPos) = double(data.images(i,j,pos(actresPos)))/double((pAct(2) + (pos(actresPos)-(Options.start_stim + Options.dur_stim*((nb-1)*2))+1)*dt *pAct(1)));
                                    posAct(actresPos) = (double(posAct(actresPos)) * 100) - 100;
                                %if rest
                                elseif ((pos(actresPos) >= (Options.start_stim + Options.dur_stim*(((nb-1)*2)+1)+1)) && (pos(actresPos) <=  (Options.start_stim + Options.dur_stim*(nb)*2)))
                                    posRest(actresPos) = double(data.images(i,j,pos(actresPos)))/double((pRest(2) + ((pos(actresPos)-((Options.start_stim + Options.dur_stim*(((nb-1)*2)+1))+1))*dt *pRest(1))));
                                    posRest(actresPos) = (double(posRest(actresPos)) * 100) - 100;
                                end
                            end
                            %squeeze the arrays and assing them to the act
                            %and rest position of the array.
                            posAct(posAct==0) = [];
                            posRest(posRest==0) = [];
                            activation.meanValue(i,j) = mean(posAct);
                            activation.stdValue(i,j) = std(posAct);
                            deactivation.meanValue(i,j) = mean(posRest);
                            deactivation.stdValue(i,j) = std(posRest);
                            %deactivation
                        end
                    end
                    clear signal_interval;
                    clear pos;
                    clear mean_signal2;
                    clear mean_signal;
                    clear peak_total;
                end
            end
        end

    %Signal case
    case ('Signal')
        N = size(data.signal,2);% Number of files
        dF = Fs/N;% Frequency spacing
        f = -Fs/2:dF:Fs/2-dF;% Frequency scale  
        BPF = ((0.6 < abs(f)) & (abs(f) < 1.5));
        
        signal_cardiac(1,:) =  double(data.signal);
        signal_cardiacaux(1,:)=detrend(signal_cardiac,'constant');      
        spektrum = ifftshift (abs(fft(signal_cardiacaux)));
        spektrum = BPF.*spektrum;
        [x,y] = find(spektrum==max(spektrum));         
        cardiac_freq = dF*(size(f,2)/2-y(1,1)+1);
        OptionsPeak.block_size = round((1/cardiac_freq)/dt);
        blocks = round(size(signal_cardiac,2)/OptionsPeak.block_size);
        OptionsPeak.ini = 0;
        OptionsPeak.estimatedInterval = round(300/(dt*1000))+1;%Variations up to 200ms in time
        mean_signal = double(mean(signal_cardiacaux(signal_cardiacaux>0)));
        signal_interval = signal_cardiacaux((1):(1+OptionsPeak.block_size));
        %Location of the first peak
        [pos(1),~,ini] = peakdet(signal_interval, (1.5*mean_signal),OptionsPeak);                 
        if ini <= OptionsPeak.estimatedInterval
        	OptionsPeak.ini = ini+OptionsPeak.block_size;
            pos(1) = 0;
        else
            OptionsPeak.ini = ini;
        end
        clear signal_interval;
        %Search for the rest of the peaks
        for b = 1:blocks-10
        	signal_interval = signal_cardiac((OptionsPeak.ini-OptionsPeak.estimatedInterval):(OptionsPeak.ini+OptionsPeak.estimatedInterval));
            [pos(b+1),~,OptionsPeak.ini] = peakdet(signal_interval,(1.5*mean_signal),OptionsPeak);
        end
        value_total = 0;
        pos(pos==0) = [];
        %Calculate the total mean peak
        mean_signal2 = mean(signal_cardiac(signal_cardiac>0));
        for k = 1:size(pos,2)        
        	value_total = double(value_total) + double(signal_cardiac(1,pos(k)));
            peak_percentage(k) = ((signal_cardiac(1,pos(k)) - mean_signal2)/mean_signal2)*100;
        end        
        mean_peak = double(value_total)/size(pos,2);
        peak_total = (double(mean_signal2)*(size(signal_cardiac,2)/(N-size(pos,2))))- (double(value_total)/(N-size(pos,2)));
        %peak_height = double(mean_peak) - double(peak_total);
        peak_height = mean(peak_percentage);
        std_value = std(peak_percentage);
        %std_value = std(signal_cardiac(pos));
        peak_positions = pos;
        clear signal_interval;
        clear pos;
        clear mean_signal2;
        clear mean_signal;
        clear peak_total;
    case '2DPeaks'
        
        peak_height = zeros(size(mask));%Initialize mean array
        std_value = zeros(size(mask));%Initialize std array
        pos = Options.positions;
        if isfield(Options,'increase')
            increase = Options.increase;
        else
            increase = 0;
        end
        for i=1:size(data.images,1)
            for j=1:size(data.images,2)
                if mask(i,j)~=0
                    signal_cardiac(1,:) =  double(data.images(i,j,:));
                    N = size(signal_cardiac,2);
                    value_total = 0;
                    %Calculate the total mean peak
                    for k = 1:size(pos,2)        
                        value_total = double(value_total) + signal_cardiac(1,(pos(k)+increase));
                    end
                    mean_signal2 = double(mean(signal_cardiac(signal_cardiac>0)));%Mean value of the whole signal
                    mean_peak = double(value_total)/size(pos,2);%Mean value for all the peaks

                    approximateMean = (double(mean_signal2)*(size(signal_cardiac,2)/(N-size(pos,2))))- (double(value_total)/(N-size(pos,2)));%Aprox. mean without peaks
                    peak_height(i,j) = double(mean_peak) - double(approximateMean);
                    %peak_height(i,j) = double(mean_peak) - min(signal_cardiac);
                    std_value(i,j) = std(signal_cardiac(pos));%Std of the peaks
                    peak_positions(i,j).pos = pos; 
                    clear signal_cardiac;
                end
            end
        end
end
