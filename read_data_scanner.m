function data = read_data_scanner(Options)
%
%
%Function which is going to return the data for every scanner, and which is
%going to be processed afterwards
%
%data will be and structured containing the parameters (data.parameters),
%the images (data.images) and the type(data.type), and if it is required
%the std and the mean value of the images. (data.std, data.mean)
%
%
% V 1.0 Creation of the file 30.11.2011 by David Lopez
% V 1.1 If filename is given it will be move foward by David Lopez 14.12.2011
% V 1.2.Bug fix. The protocol value is properly set in case of bruker
% data 15.12.2011
% V 1.3 Bug Fix with the value of the mean_std paramater in case of Bruker
% data 16.12.2011
% V 1.4 Bug Fix Now the value of the Mean_std variable expects True or
% False value David Lopez 19.12.2011
% V 1.5 Bug Fix the number of slices is given to calculate the mean and std
% values by David Lopez 16/01/2012
% V 1.6 Multishot case also added by David Lopez 23.03.2012
% V 1.7 Dicom data is also accepted by David Lopez 28.05.2012
% V 1.8 Bug Fix. in Bruker Data by David Lopez 12.06.2012
% V 1.9 The order parameter has been introduced to indicate if the data
% slices are interleaved or noninterleaved by David Lopez 13.06.2012
% V 2.0 Now it accepts the parameters sk_ku which is going to determined if
% the skewness and kurtosis of the data has to be calculated when the std
% is performed by David Lopez 02.07.2012
% V 2.1 The dicom parameters have been updated by David Lopez 14.07.2012
% V 2.2 A parameter to discard initial slices  in the std calculations has been introduced
%  by David Lopez 17.08.2012
% V 2.3 The crop slices has been uptated to crop the images also to the
% initial data load by David Lopez 23.08.2012
% V 2.4 Bug Fix. Now the crop images is working in Nifti files 10.10.2012
% by David Lopez
% V 2.5 Bug. Fix. When we receive NIFTI files the parameters to calculate
% the std and mean are read now correctly by David Lopez 31.1.2013
% V 2.6 Bug Fix. When the Nifti files arrive we have introduce the
% load_untouch_nii m,ethod to recover the data so the images = images.img
% from the method by David Lopez 11.04.0213
% V 2.7 The mask can be given now in the Options struct to avoid
% calculations by David Lopez 29.04.2013
% V 2.8 Bug Fix. The dicom data has to be reslice because all the slices
% are loaded one after the other by David Lopez 15.05.2013
% V 2.9 Adaptation to the new Dicom data where the NewDicom files are read
% in the same way as the Dicom ones by David Lopez 29.07.2015


if nargin == 0
    type = 'Bruker';
end

if nargin == 1 && isempty(Options.type)
    type = 'Bruker';
else
    type = Options.type;
    if strcmp(type,'NewDicom')
        type = 'Dicom';
    end
end

if isfield(Options,'filename')
    filename = Options.filename;
else
    filename = '';
end

if nargin == 1 && ~isfield(Options,'order')
    order = 'Interleaved';
    if isfield(Options,'shot')
        shot = Options.shot;
    else
        shot = 'shot';
    end
else   
    order = Options.order;
    if strcmp(order,'Interleaved')
        shot = 'shot';
    else
        shot = 'Multishot';
    end
end

if nargin == 1 && ~isfield(Options,'mean_std')
    mean_std = 'true';   
else
    mean_std = Options.mean_std;
end

if nargin == 1 && ~isfield(Options,'protocol') && strcmp(type,'Philips')
    protocol = 'Nifti';
elseif strcmp(type,'Bruker') || strcmp(type,'Dicom')
    protocol = '';
else
    protocol = Options.protocol;
end

if nargin == 1 && ~isfield(Options,'sk_ku')
    sk_ku = 'yes';   
else
    sk_ku = Options.sk_ku;
end

if nargin == 1 && ~isfield(Options,'crop_slices')
    crop_slices = 0 ;
else
    crop_slices = Options.crop_slices;
end

switch type
    case 'Bruker'
        fprintf('Reading Bruker data \n');
        %Get the file path
        if nargin == 1 && isempty(Options.filename)
            [fname,pname] = uigetfile('*.*', 'Select the 2dseq data');
            if pname
                filename = strcat(pname,fname);
            end
            %Set the path folder
            cd(pname);               
        end

        %It Returns the images, parameters and type in case of Bruker data
        images = process_images(filename);
        %It Returns the mean value and the std value of the data
        if strcmp(mean_std,'True') || strcmp(mean_std,'true') || strcmp(mean_std,'Yes') || strcmp(mean_std,'yes')
            Optionsmean.type = type;
            if  ~isfield(Options,'scans')
                scans = images.parameters.NR;
            else
                scans = Options.scans;
            end
            if  ~isfield(Options,'echoes')
                echoes = images.parameters.NI;
            else
                echoes = Options.echoes;
            end
            Optionsmean.scans = scans;
            Optionsmean.echoes = echoes;
            Optionsmean.slices = images.parameters.ns;
            Optionsmean.order  = order;
            fprintf('Calculating the standard deviation and mean value \n');
            error = validate_dimensions(scans, echoes,Optionsmean.slices,size(images.data,3)*size(images.data,4));
            if error == 1
                errordlg('The number of echoes or the number of scans is incorrect. It doesnt correspond with the read file');
                return;
            end
            Optionsmean.sk_ku = sk_ku;
            Optionsmean.crop_slices = crop_slices;
            if strcmp(sk_ku,'Yes') || strcmp (sk_ku,'yes')
                [mean_value,std_value,skewness,kurtosis] = calculate_std_mean(images.data,Optionsmean);
                data.kurtosis = kurtosis;
                data.skewness = skewness;
                data.mean = mean_value;
                data.std = std_value;   
            else
                [mean_value,std_value] = calculate_std_mean(images.data,Optionsmean);
                data.mean = mean_value;
                data.std = std_value;   
            end
        end
        
        data.type = images.type_data;
        data.images = images.data;
        data.parameters = images.parameters;
        data.error = images.error;

    case 'Philips'
        fprintf('Reading Philips data');
        if strcmp(protocol,'PAR')
            %header = read_philips_header();
            OptionsRead.protocol = protocol;
            OptionsRead.filename = filename;
            [folder,filen]=fileparts(filename);
            filenameREC=fullfile(folder,[filen '.rec']);
            OptionsRead.filenameREC = filenameREC;
            OptionsRead.shot = shot;
            [images,parameters] = read_philips_data(OptionsRead);
        elseif strcmp(protocol,'Nifti')
            OptionsRead.protocol = protocol;
            OptionsRead.filename = filename;
            [imagesaux,parameters] = read_philips_data(OptionsRead);
            if isfield(imagesaux,'img') 
                images = imagesaux.img;
            else
                images = imagesaux;
            end
        else
            errordlg('This protocol is not supported');
        end
        %Call to the calculation of std and mean values if option is set to Yes
        if strcmp(mean_std,'True') || strcmp(mean_std,'true') || strcmp(mean_std,'Yes') || strcmp(mean_std,'yes')
            Optionsmean.type = type;
            if strcmp(protocol, 'PAR')
                Optionsmean.slices = str2num(parameters.MaxnumberofslicesLocations);
                Optionsmean.scans = str2num(parameters.Maxnumberofdynamics);
                Optionsmean.echoes = str2num(parameters.Maxnumberofechoes);
            else
                Optionsmean.slices = parameters.Dimensions(3);
                Optionsmean.scans = parameters.Dimensions(4);
                Optionsmean.echoes = parameters.Dimensions(5);
            end
            if isfield(Options,'Mask')
                Optionsmean.Mask = Options.Mask;
            end
            Optionsmean.protocol = protocol;
            Optionsmean.order = order;
            fprintf('Calculating the standard deviation and mean value \n');
            %error = validate_dimensions(Options.scans, Options.echoes,Optionsmean.slices,size(images,3)*size(images,4));
            %if error == 1
            %    errordlg('The number of echoes or the number of scans is incorrect. It doesnt correspond with the read file');
            %    return;
            %end            
            Optionsmean.sk_ku = sk_ku;
            Optionsmean.crop_slices = crop_slices;
            if strcmp(sk_ku,'No') ||strcmp(sk_ku,'no')
                [mean_value,std_value] = calculate_std_mean(images,Optionsmean);
                data.mean = mean_value;
                data.std = std_value;
            else
                [mean_value,std_value,skew,kurt] = calculate_std_mean(images,Optionsmean);
                data.mean = mean_value;
                data.std = std_value;
                data.skew = skew;
                data.kurt = kurt;
            end
        end
        
        data.type = type;
        data.images = images;
        data.parameters = parameters;

    case 'Dicom'
        %Get the file path if the user hasnt introduced it
        if nargin == 1 && isempty(Options.filename)
            [fname,pname] = uigetfile('*.*', 'Select the dicom file');
            if pname
                filename = strcat(pname,fname);
            end
        end
        
        data.images = dicomread(filename);
        data.parameters = dicominfo(filename);
        %for i = 1:data.parameters.Private_2001_1018
        %    aux(:,:,i,:) = data.images(:,:,1,(((i-1)*(data.parameters.NumberOfFrames/data.parameters.Private_2001_1018))+1):((data.parameters.NumberOfFrames/data.parameters.Private_2001_1018)*i));
        %end
        %data = rmfield(data,'images');
        %data.images = aux;
        %{
        if strcmp(mean_std,'True') || strcmp(mean_std,'true') || strcmp(mean_std,'Yes') || strcmp(mean_std,'yes')
            fprintf('Reading dicom data');
            Optionsmean.type = type;
            if  ~isfield(Options,'scans')
                Optionsmean.scans = data.parameters.NumberOfFrames /(data.parameters.Private_2001_1014*data.parameters.Private_2001_1018);
            else
                Optionsmean.scans = Options.scans;
            end
            if ~isfield(Options,'echoes')
                Optionsmean.echoes = data.parameters.Private_2001_1014;
            else
                Optionsmean.echoes = Options.echoes;
            end
            Optionsmean.slices = data.parameters.Private_2001_1018;
            Optionsmean.order = order;
            error = validate_dimensions(Optionsmean.scans, Optionsmean.echoes, Optionsmean.slices,size(data.images,3)*size(data.images,4));
            if error == 1
                    errordlg('The number of echoes or the number of scans is incorrect. It doesnt correspond with the read file');
                    return;
            end    
            Optionsmean.crop_slices = crop_slices;
            Optionsmean.sk_ku = sk_ku;
            if strcmp(mean_std,'True') || strcmp(mean_std,'true') || strcmp(mean_std,'Yes') || strcmp(mean_std,'yes')
                if strcmp(sk_ku,'No') ||strcmp(sk_ku,'no')
                    fprintf('Calculating the standard deviation and mean value \n');
                    [mean_value,std_value] = calculate_std_mean(data.images,Optionsmean);
                    data.mean = mean_value;
                    data.std = std_value;
                else
                    fprintf('Calculating the standard deviation, mean value, skewness and kurtosis \n');
                    if strcmp(shot,'shot')
                        aux = zeros([size(data.images,1) size(data.images,2) Optionsmean.slices size(data.images,4)/Optionsmean.slices]);
                        for i=1:size(aux,3)
                            for j=1:size(aux,4)
                                 aux(:,:,i,j) = data.images(:,:,j + size(aux,4)*(i-1));
                            end
                        end
                        data.images = aux;
                    else
                        data.images = reshape(data.images,[size(data.images,1) size(data.images,2) Optionsmean.slices size(data.images,4)/Optionsmean.slices]);
                    end
                    [mean_value,std_value,skewness,kurtosis] = calculate_std_mean(data.images,Optionsmean);
                    data.mean = mean_value;
                    data.std = std_value;
                    data.skewness = skewness;
                end
            end
        end
        %}
        
end

%Crop the initial slices if it has been specified
if crop_slices > 1
	dataaux = data.images;
    clear data.images;
    if ndims(dataaux) == 3 
        data.images = dataaux(:,:,crop_slices+1:end);
    elseif ndims(dataaux) == 4 
        %if  ~strcmp(protocol,'Nifti')
            data.images = dataaux(:,:,:,crop_slices+1:end);
        %else
         %   data.images.img = dataaux.img(:,:,:,crop_slices+1:end);
        %end
    end
end

end

function error = validate_dimensions(scan, echoes,slices,total)
    calculation = scan*echoes*slices;
    if total~=calculation;
        error = 1;
    else
        error = 0;
    end
end