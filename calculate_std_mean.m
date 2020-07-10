function [mean_value,std_value,skew,kurt]=calculate_std_mean(array,Options)

% Method which is going to return the standard deviation and the mean value
% of an array of values
%
% [A,B] = calculate_std_mean(C) It will return the mean value A and the
% standard deviation B from array of valies C.
% [A,B] calculate_std_mean(C,Options)--> Same as before but now the options
% are given just in case Philips data wants to be analyzed
%
% V 1.0 Creation of the document -> David Lopez
% V 1.1 OPtion type included for the possibility of using Philip data by David Lopez 1.12.2011
% V 1.2 Adapted to multislice by David Lopez 17/1/2012
% V 1.3 Bug fix --> The value of Options. slices, scans and echoes is checked at the
% beginning  by David Lopez 02/02/2012
% V 1.4 Bug fix --> Case Philips2 has been introduced to treat the case in
% which the data is not given echo1-echo2-echo3-echo4-echo5-echo6 but
% echo1-echo1-echo2-echo2.. by David Lopez 07/02/2012
% V 1.5 Bug fix --> The default value of scans has been corrected by David
% Lopez 07/02/2012
% V 1.6 Bruker processing removed and names of the type of processing
% changed from Philips and Philips yo Interleaved and Noninerleaved by
% David Lopez 13.06.2012
% V 1.7 The skewness and the kurtosis of the data is calculated to
% determine the "Gaussianity" of the data  by David Lopez 02.07.2012
% V 1.8 Now we calculate the mask to perform a cropping on the imput data
% and reduce the calculation time by David Lopez 12.07.2012
% V 1.9 Introduced the parameters to eliminate the first slices of a
% sequence by David Lopez 17.08.2012
% V 2.0 Parallel programming has been introduced to improve the processing
% speed by David Lopez 14.09.2012
% V 2.1 Bug fix in the parallel programming when the number of dimensions
% was equal to 3 in the single slice case by David Lopez 04.12.2012
% V 2.2 Bug fix in the 4D case because in the case of PAR REC data the
% slice number comes in the third position while in Nifti it does in the
% fourth one by David Lopez 08.02.2013
% V 2.3 Now the mask can be given in the Options file by David Lopez
% 29.04.2013

if nargin==0
    errordlg('There is no sufficient input values ');    
end

if nargin < 2
    order = 'Interleaved';
else
    order = Options.order;
end

if nargin < 2 || ~isfield(Options,'slices')
    slices = size(array,3);
else
    slices = Options.slices;
end

if nargin < 2 || ~isfield(Options,'echoes')
    echoes = size(array,4)/100;
else
    echoes = Options.echoes;
end

if nargin < 2 || ~isfield(Options,'scans')
    scans = size(array,4)/echoes;
else
    scans = Options.scans;
end

if nargin < 2 || ~isfield(Options,'sk_ku')
    sk_ku = 'No';
else
    sk_ku = Options.sk_ku;
end

if nargin < 2 || ~isfield(Options,'crop_slices')
    crop_slices = 0 ;
else
    crop_slices = Options.crop_slices;
    scans = scans-crop_slices;
end

arrayMask = [];

if size(array,4) > 1
    if ~isfield(Options,'Mask')
        array_aux2 = array(:,:,:,(crop_slices + 1):end);
        clear array;
        array = array_aux2;
        OptionsMask.limit = 100;
        OptionsMask.neighbors = 8;
        OptionsMask.neighborstype2 = 8;
        OptionsMask.dotype2correction = true;       
        mask = calculate_mask(array(:,:,:,1),OptionsMask);
        totalMask = ~mask;
    else 
        totalMask = Options.Mask;
    end
    
    [totalMask, arrayMask] = zelaouscrop(totalMask);
    arrayMask(1,1) = arrayMask(1,1);
    arrayMask(2,1) = arrayMask(2,1);
    arrayMask(1,2) = arrayMask(1,2);
    arrayMask(2,2) = arrayMask(2,2);
    %if size(array,3)>1
    %    arrayMask(3,1) = arrayMask(3,1)-2;
    %    arrayMask(3,2) = arrayMask(3,2)+2;
    %end
    for i=1:size(array,4)
        aux_vector(:,:,:,1) =  array(:,:,:,i);
        array_aux(:,:,:,i) = zelaouscrop(aux_vector,arrayMask);
    end
else
    if ~isfield(Options,'Mask')
        array_aux2 = array(:,:,(crop_slices + 1):end);
        clear array;
        array = array_aux2;
        OptionsMask.limit = 20;
        OptionsMask.neighbors = 8;
        OptionsMask.neighborstype2 = 8;
        OptionsMask.dotype2correction = true;
        mask = calculate_mask(array(:,:,1),OptionsMask);
        totalMask = ~mask;
    else
        totalMask = Options.Mask;
    end
    [totalMask, arrayMask] = zelaouscrop(totalMask);
    arrayMask(1,1) = arrayMask(1,1)-1;
    arrayMask(2,1) = arrayMask(2,1)-1;
    arrayMask(1,2) = arrayMask(1,2)+1;
    arrayMask(2,2) = arrayMask(2,2)+1;
    for i=1:size(array,3)
        aux_vector =  array(:,:,i);
        array_aux(:,:,i) = zelaouscrop(aux_vector,arrayMask);
    end
end
clear array;
array = array_aux;

%Calculation of mean and std
switch order
    case 'Interleaved'
    %Process data for 1 slice
    if(slices == 1)
        parfor x=1:size(array,1)
            aux(x,:,:,:,:) = Interleaved1array(array(x,:,:,:,:),scans,echoes);
        end
        parfor x=1:size(aux,1)            
            if strcmp(sk_ku,'Yes')||strcmp(sk_ku,'yes')
                [mean_value(x,:,:,:,:),std_value(x,:,:,:,:),skew(x,:,:,:,:),kurt(x,:,:,:,:)] = Interleaved1calc(aux(x,:,:,:,:),sk_ku);
            else
                [mean_value(x,:,:,:,:),std_value(x,:,:,:,:)] = Interleaved1calc(aux(x,:,:,:,:),sk_ku);
            end
        end
    %Process the data for several slices
    else
        parfor x=1:size(array,1)
            aux(x,:,:,:,:) = Interleaved1mult_array(array(x,:,:,:),slices,scans,echoes);
        end
        parfor x=1:size(aux,1)
            if strcmp(sk_ku,'Yes')||strcmp(sk_ku,'yes')
                [mean_value(x,:,:,:,:),std_value(x,:,:,:,:),skew(x,:,:,:,:),kurt(x,:,:,:,:)] = Interleaved1multcalc_array(aux(x,:,:,:,:),sk_ku);
            else
                [mean_value(x,:,:,:,:),std_value(x,:,:,:,:)] = Interleaved1calc(aux(x,:,:,:,:),sk_ku);
            end
        end
    end
    %This case is created in case that data comes in another order. First
    %case for multislice.
    case 'Noninterleaved'
        parfor x=1:size(array,1)
            aux(x,:,:,:,:) = NonInterleaved1mult_array(array(x,:,:,:,:),slices,scans,echoes);
        end
        parfor x=1:size(aux,1)
            if strcmp(sk_ku,'Yes')||strcmp(sk_ku,'yes')
                [mean_value(x,:,:,:,:),std_value(x,:,:,:,:),skew(x,:,:,:,:),kurt(x,:,:,:,:)] = NonInterleaved1multcalc_array(aux(x,:,:,:,:),sk_ku);
            else
                [mean_value(x,:,:,:,:),std_value(x,:,:,:,:)] = Interleaved1calc(aux(x,:,:,:,:),sk_ku);
            end
        end
    otherwise
        errordlg();
end
if size(arrayMask,1) > 2
    for i = 1: size(mean_value,4)
        aux_mean(:,:,:,i) = zelaouscrop(mean_value(:,:,:,i),arrayMask,[],true);
        aux_std(:,:,:,i) = zelaouscrop(std_value(:,:,:,i),arrayMask,[],true);   
    end
else
    for i = 1: size(mean_value,3)
        aux_mean = zelaouscrop(mean_value,arrayMask,[],true);
        aux_std = zelaouscrop(std_value,arrayMask,[],true);   
    end
end      
clear mean_value;
clear std_value;
mean_value = aux_mean; 
std_value = aux_std;
if strcmp(sk_ku,'Yes')||strcmp(sk_ku,'yes')
    if size(arrayMask,1) > 2
        for i = 1: size(skew,4)
            skew_aux(:,:,:,i) = zelaouscrop(skew(:,:,:,i),arrayMask,[],true);
            kurt_aux(:,:,:,i) = zelaouscrop(kurt(:,:,:,i),arrayMask,[],true);
        end
    else
        for i = 1: size(skew,4)
            skew_aux(:,:,i) = zelaouscrop(skew(:,:,i),arrayMask,[],true);
            kurt_aux(:,:,i) = zelaouscrop(kurt(:,:,i),arrayMask,[],true);
        end
    end
    clear skew;
    clear kurt;
    skew = skew_aux;
    kurt = kurt_aux;
end
end

function aux = Interleaved1array(array,scans,echoes)
	for y=1:size(array,2)
    	value = 0;
        for z=1:scans
        	for w=1:echoes
            	value = value + (w-(w-1));
                aux(y,z,w) = array(1,y,value);
            end
        end
	end
end

function [mean_value,std_value,skew,kurt] = Interleaved1calc(aux,sk_ku)
    if ndims(aux) == 4
        for y=1:size(aux,2)
            for i=1:size(aux,3)
                mean_value(y,i,1)=mean(aux(1,y,i,:));
                std_value(y,i,1)=std(double(aux(1,y,i,:)),1);
                if strcmp(sk_ku,'Yes')||strcmp(sk_ku,'yes')
                    skew(y,i,1) = skewness(double(aux(1,y,i,:)));
                    kurt(y,i,1) = kurtosis(double(aux(1,y,i,:)));                        
                end
            end
        end
    else
        for y=1:size(aux,2)
            mean_value(y,1)=mean(aux(1,y,:));
            std_value(y,1)=std(double(aux(1,y,:)),1);
            if strcmp(sk_ku,'Yes')||strcmp(sk_ku,'yes')
                skew(y,1) = skewness(double(aux(1,y,:)));
            	kurt(y,1) = kurtosis(double(aux(1,y,:)));                        
            end
        end
    end
end
        
function aux = Interleaved1mult_array(array,slices,scans,echoes)
for y=1:size(array,2)                
	for r=1:slices
    	value = 0;
        	for z=1:scans
            	for w=1:echoes
                	value = value + (w-(w-1));
                    aux(y,r,z,w) = array(1,y,r,value);
                end
            end
    end
end
end

function [mean_value,std_value,skew,kurt] = Interleaved1multcalc_array(aux,sk_ku)    
    for y=1:size(aux,2)
    	for r=1:size(aux,3)
        	for i=1:size(aux,5)
            	mean_value(y,r,i,1)=mean(aux(1,y,r,i,:));
                std_value(y,r,i,1)=std(double(aux(1,y,r,i,:)));
                if strcmp(sk_ku,'Yes')||strcmp(sk_ku,'yes')
                    skew(y,r,i,1) = skewness(double(aux(1,y,i,:)));
                    kurt(y,r,i,1) = kurtosis(double(aux(1,y,i,:)));
                end
            end
        end
    end
end

function aux = NonInterleaved1mult_array(array,slices,scans,echoes)
for y=1:size(array,2)                
	for r=1:slices
    	for w=1:echoes
        	for z=1:scans
            	aux(y,r,z,w) = array(1,y,r,z);%Need to be modified for multiecho
            end
        end
    end
end
end

function [mean_value,std_value,skew,kurt] = NonInterleaved1multcalc_array(aux,sk_ku)    
    for y=1:size(aux,2)
        for r=1:size(aux,3)
            %for i=1:size(aux,4)
                mean_value(y,r,1)=mean(aux(1,y,r,:));
                std_value(y,r,1)=std(double(aux(1,y,r,:)));
                if strcmp(sk_ku,'Yes')||strcmp(sk_ku,'yes')
                	skew(y,r,1) = skewness(double(aux(1,y,r,:)));
                    kurt(y,r,1) = kurtosis(double(aux(1,y,r,:)));
                end
        	%end
    	end
    end
end
