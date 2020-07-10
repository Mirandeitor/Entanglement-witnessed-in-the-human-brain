function [mask] = calculate_mask(data_array,Options)

%Function which is going to calculate the approximate data of the brain to
%avoid too long calculations and eliminate data outside the brain. (Based on the CTM
%method but using in this case only the magnitude data).

% [mask] = calculate_mask(A,Options) It's gonna calculate the mask of the
% data array A using some established options.

% V 1.0 Creation of the document by David Lopez 05/02/2012 
% V 1.1 Adaptation of the method to Bruker data by David Lopez 28/02/2012
% V 1.2 A parfor loop has been added to increase speed by David Lopez
% 14.09.2012
% V 1.3 A snr_mask has been added to improve the determination of the brain
% mask by David Lopez 20.12.2012

%Check the input parameters
%fprintf('Calculating mask \n');

AUX_COUNTER = 0;

if nargin == 0 
    error('No input data has been provided');
end

if nargin == 1 
    mag_limit = 15;
    if size(dataArray,3) > 1
        mag_neighbors = 16;
        mag_neighbors_type2 = 17;
    else
        mag_neighbors = 7;
        mag_neighbors_type2 = 9;
    end
else
    if ~isfield(Options,'limit')
        mag_limit = 15;
    else
        mag_limit = Options.limit;
    end
    if ~isfield(Options,'neighbors')
        mag_neighbors = 16;
    else
    mag_neighbors = Options.neighbors;
    end
    if ~isfield(Options,'neighborstype2')
        mag_neighbors_type2 = 7;
    else
        mag_neighbors_type2 =  Options.neighborstype2;
    end
end

%Perform the mask calculation
magnitudeMask = data_array;
magnitudeMask(magnitudeMask <= mag_limit) = false;
magnitudeMask(magnitudeMask > mag_limit) = true;

%Reduce the number of errors analyzing the neighbors
if mag_neighbors ~= 0
    %fprintf('Reducing type 1 error in magnitude data \n');
    finalmask = checkneighbors(~magnitudeMask,double(magnitudeMask),mag_neighbors,AUX_COUNTER);
end

%Reduce the type II errors
if Options.dotype2correction
     %fprintf('Reducing type 2 error in magnitude data \n');
     finalmask = checkneighbors(~finalmask,double(finalmask),mag_neighbors_type2,AUX_COUNTER);
end

mask = ~finalmask;
%Calculate a snr_mask to perform a more accurate masking

automaticSnrRoiParameter = false;
snrMin  = 0.1;
homoRoiBoolean = false;
roi = finalmask;
if isempty(roi)
    if ndims(data_array) == 3
        roi(1:round(end*automaticSnrRoiParameter)+2,1:round(end*automaticSnrRoiParameter)+2,1:round(end*automaticSnrRoiParameter)+2) = true; % +2 for the case that round(end*automaticSnrRoiParameter) is zero
    elseif ndims(data_array) == 2
        roi(1:round(end*automaticSnrRoiParameter)+2,1:round(end*automaticSnrRoiParameter)+2) = true; % +2 for the case that round(end*automaticSnrRoiParameter) is zero
    else
        error('This data matrix size is currently not supported (only 2D and 3D)')
    end
end
if ~homoRoiBoolean
    noiseLevel = mean(data_array(roi));
else
    noiseLevel = std(data_array(roi));
end

threshold = snrMin * noiseLevel;

mask_snr = data_array >= threshold;

fprintf('End of calculating mask \n');
mask = ~(finalmask .* mask_snr);
%}

end

%Function which is going to analyze the neighbours
function [finalmask] = checkneighbors(inputmask,inputdatamask,numberofconnectivity,AUX_COUNTER)
%This function compute the connectivity criteria for the input data
%inputmask = Input data is the data to be checked
%inputdatamask = Input thresholded mask
%numberofconnetivity = value which indicates the number of neighbors that
%have to be reliable to restore the pixel.

szindata = size(inputmask);
numberofDimensions = ndims(inputmask);

%Create the shift vector to identify the surrounding neighbors
switch numberofDimensions
    case 2
        numberofbeighborhoodvoxels = 8;
        shiftVector(:,1) = [1 0];
        shiftVector(:,2) = [-1 0];
        shiftVector(:,3) = [0 1];
        shiftVector(:,4) = [0 -1];
        shiftVector(:,5) = [1 1];
        shiftVector(:,6) = [-1 1];
        shiftVector(:,7) = [1 -1];
        shiftVector(:,8) = [-1 -1];
    case 3
        numberofbeighborhoodvoxels = 26;
        shiftVector = zeros(numberofDimensions, numberofbeighborhoodvoxels);
        % direct neighbors
        shiftVector(:,1)  = [ 1,  0,  0];
        shiftVector(:,2)  = [-1,  0,  0];
        shiftVector(:,3)  = [ 0,  1,  0];
        shiftVector(:,4)  = [ 0, -1,  0];
        shiftVector(:,5)  = [ 0,  0,  1];
        shiftVector(:,6)  = [ 0,  0, -1];
                
        % neighbors right to x0, y0, z0
        shiftVector(:,7)  = [ 1,  1,  0];
        shiftVector(:,8)  = [ 1,  1,  1];
        shiftVector(:,9)  = [ 1,  1, -1];
        shiftVector(:,10) = [ 1, -1,  0];
        shiftVector(:,11) = [ 1, -1,  1];
        shiftVector(:,12) = [ 1, -1, -1];
        shiftVector(:,13) = [ 1,  0,  1];
        shiftVector(:,14) = [ 1,  0, -1];
                
        % neighbors left to x0, y0, z0
        shiftVector(:,15) = [-1,  1,  0];
        shiftVector(:,16) = [-1,  1,  1];
        shiftVector(:,17) = [-1,  1, -1];
        shiftVector(:,18) = [-1, -1,  0];
        shiftVector(:,19) = [-1, -1,  1];
        shiftVector(:,20) = [-1, -1, -1];
        shiftVector(:,21) = [-1,  0,  1];
        shiftVector(:,22) = [-1,  0, -1];
                
        % neighbors up and below to x0, y0, z0
        shiftVector(:,23) = [0,  1,  1];
        shiftVector(:,24) = [0,  1, -1];
        shiftVector(:,25) = [0, -1,  1];
        shiftVector(:,26) = [0, -1, -1];
                
    otherwise
        error('The following dimension is not allowed');
end

 % Correction of errors: If we are correcting the typ of error 1, it is necessary to analyze the pixels that have been
 % removed by the thresshold method. Those whose connectivity is bigger or equal than X will be restored. if we are correcting
 % error 2 we have to analyze all values that havent been removed by the thresholding. Those whose connectivity is lower than X will be removed 
 
    connectivityMask = zeros(szindata);
    parfor i = 1:numberofbeighborhoodvoxels
        connectivityMask = connectivityMask + circshift(inputdatamask, shiftVector(:,i));
    end

    if AUX_COUNTER < 1
        inputmask(inputmask) = ~(connectivityMask(inputmask) >= numberofconnectivity);
    end
        
    if AUX_COUNTER == 1
        inputmask(~inputmask) = (connectivityMask(~inputmask) < numberofconnectivity);
    end
        
    AUX_COUNTER = AUX_COUNTER + 1;
    finalmask = ~inputmask;
end