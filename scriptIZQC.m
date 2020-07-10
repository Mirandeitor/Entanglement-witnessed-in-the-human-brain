clear all;

%Depending of the type of data and experiment the Options.type nad the Options.protocol needs to be
%changed.

Options.type= 'NewDicom';
%Options.protocol = 'PAR';
Options.mean_std = 'Yes';
Options.crop_slices = 100;
Options.filename ='/Users/lopezd/Documents/PhD Data/Data/2015/Method paper Data/Hypo Hyper ventilation test/IM_0081.dcm';
data = read_data_scanner(Options);
data.images = squeeze(data.images);
%data.parameters.Private_2001_1023
%data.images = ((double(data.images)*double(data.parameters.RescaleSlope))+data.parameters.RescaleIntercept)./double(data.parameters.RescaleSlope*data.parameters.Private_2005_100e);
%data1 = data;
%Options.filename = '/Volumes/David/Local Disk/Phd/Phd/data/2013/data from study/Young/25.06.2013/429/ab250613ek429_10_1.PAR';
%data1 = read_data_scanner(Options);

Options.type = 'Dicom';
Options.filename ='/Users/lopezd/Documents/PhD Data/Data/2015/Method paper Data/David Diaz/Frequency/6.1/DICOM/IM_0001.dcm';
[file2,name2,ext2] = fileparts(Options.filename);
for i=101:1000
    if i < 10
        Options.filename = strcat(file2,'/','IM_000',num2str(i));
    elseif i == 9 || i < 100
        Options.filename = strcat(file2,'/','IM_00',num2str(i));
    elseif i == 100 || i < 1000
        Options.filename = strcat(file2,'/','IM_0',num2str(i));
    else
        Options.filename = strcat(file2,'/','IM_',num2str(i));
    end
    dataaux = read_data_scanner(Options);
    data(:,:,(i-100)) = dataaux.images;
end

data = ((double(data)*double(dataaux.parameters.RescaleSlope))+dataaux.parameters.RescaleIntercept)./double(dataaux.parameters.RescaleSlope*dataaux.parameters.Private_2005_100e);
%{
data1.images1 = data;
data1.parameters = dataaux.parameters;
%}
%mask_aux = load_untouch_nii('C:\Phd\Phd\data\2015\Test slab distance 2\IM_1000.nii');
%mask_aux = load_untouch_nii('C:\Phd\Phd\data\2014\data for paper\2014\new coil\Visual longer stimulation\mask_test14102014_6_1.nii')
%mask = mask_aux.img(:,:,1);
%mask = erodemask(mask,3);

%{
[root,file,extension] = fileparts(Options.filename);
filename = strcat(root,'\',file,'.nii');
nii = make_nii(data1.images1,[],[],[])
save_nii(nii,filename);


%Separate in two echoes

for i=1:size(data1.images,1)
    for j=1:size(data1.images,2)
        for k=1:size(data1.images,3)/2
            data1.images1(i,j,k) = data1.images(i,j,k*2-1);            
            data1.images2(i,j,k) = data1.images(i,j,k*2);    
        end
    end
end


%Mask
mask_aux = load_untouch_nii('C:\Phd\Phd\data\2014\data for paper\2014\new coil\Christian 22.10.2014\mask_test_22102014_32_1.nii');
mask = mask_aux.img(:,:,1);
mask = erodemask(mask,2);   
%}

%data1.images1 = data1.images;
%Calculate the total signal average of the brain
mask = calculate_magnitude_mask(data.images);
%mask = calculate_magnitude_mask(data);
mask = erodemask(mask,2); 
counter = 0;
signal_total_brain = zeros([1 size(data.images,3)]);
%signal_total_brain = zeros([1 size(data,3)]);
%signal_total_brain1 = zeros([1 size(data1.images,3)]);
%signal_total_brain = zeros([1 2100]);
for i=1:size(data.images,1)
	for j=1:size(data.images,2)
%for i=1:size(data,1)
	%for j=1:size(data,2)
        if mask(i,j)~=0
            aux(1,:) = double(data.images(i,j,:));
            %aux(1,:) = double(data(i,j,:));
            signal_total_brain = signal_total_brain + aux;
            counter = counter + 1;
            %aux1(1,:) = double(data1.images(i,j,:));
            %signal_total_brain1 = signal_total_brain1 + aux1;
        end
    end
end
%{
for i=1:size(data.images,1)
	for j=1:size(data.images,2)
        std_image(i,j) = std(double(data.images(i,j,:)));
    end
end
%}
signal_total_brain= signal_total_brain / counter;
tScale = 0.045:0.045:(0.045*6500);
stringComponent = 'Signal = {';
for i=1:size(tScale,2)
    stringComponent = strcat(stringComponent,'{',num2str(tScale(1,i)),',',num2str(signal_total_brain(1,i)),'},');
end
stringComponent(end) = '}';
stringToWrite{1} = stringComponent;

fileID = fopen('Pixel4.nb','w');
for i = 1:size(stringToWrite,2)
    fprintf(fileID,'%s\n',stringToWrite{i});
end
fclose(fileID);


%signal_cardiac(1,:) = double(data1.images(55,31,:));
%signal_cardiac = double(signal_cardiac + 300);
time_scale = 0.060:0.060:(0.060*2900);
figure;plot(time_scale,signal_cardiac,'r','linewidth',1)
hold on;
plot(double(signal_total_brain),'linewidth',1);
xlabel('Time (seconds)');
ylabel('Signal Intensity');
xlim([60 70]);
legend('Vein Signal','Tissue Signal')
print('signalVein', '-dpng', '-r600');
%Calculate peaks
dataaux.signal = signal_total_brain(1,1:end);
dataaux.parameters = data.parameters;
OptionsPeaks.DataType ='NewDicom';
dt = data.parameters.AcquisitionDuration/(100+size(signal_total_brain,2));% Time spacing
[meanHeight,stdHeight,peaksPosition] = calculatePeakStatistics(dataaux, mask, OptionsPeaks);
meanSignal = mean(signal_total_brain(1,1:end));
stdSignal = std(signal_total_brain(1,1:end));


%signal_total_brain1 = signal_total_brain1 / counter;

%{
%time_scale = (str2num(data.parameters.Repetitiontimems)/1000):(str2num(data.parameters.Repetitiontimems)/1000):size(signal_total_brain,2)*(str2num(data.parameters.Repetitiontimems)/1000);
figure;pic = plot(time_scale(1,1:900),signal_total_brain(1,1:900));
hold on;
plot(time_scale(1,1:900),signal_total_brain1(1,1:900));
xlabel('Time(seconds)');
ylabel('Signal Intensity');
%dt = str2num(data.parameters.Repetitiontimems)/1000;% Time spacing
Fs = 1/.06;% Frequency spacing
N = size(signal_total_brain{18},2);% Number of files
dF = Fs/4096;% Frequency spacing
f = -Fs/2:dF:Fs/2-dF;% Frequency scale  
BPF = ((0.6 < abs(f)) & (abs(f) < 1.5));
signal_total_brain2(1,:)=detrend(double(signal_total_brain{18}),'constant');      
spektrum = ifftshift (abs(fft(double(signal_total_brain2),4096)));
signal_total_brain(1,:)=detrend(double(signal_total_brain1),'constant');      
spektrum1 = ifftshift (abs(fft(double(signal_total_brain1),4096)));
%spektrum2 = ifftshift (abs(fft(double(my_VAR_1))));
params.Fs = Fs;
[Pxx,F] = mtspectrumc(signal_total_brain,params);
[Pxx1,F] = mtspectrumc(signal_total_brain1,params);
figure;plot(f,Pxx);
hold on
plot(f,Pxx1);
xlabel('Frequency (Hz)');
ylabel('Frequency Intensity');
xlim([0 3.5])  
%figure;plot(signal_total_brain);
%}

