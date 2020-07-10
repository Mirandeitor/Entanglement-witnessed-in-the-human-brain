function mask = calculate_magnitude_mask(data)

%Function which is estimate the magnitude mask of the input data to cover
%the whole slice

%Creation of the document by David Lopez

%Validation of the input data
if nargin < 1
    error('The input data has not been established');
end


%Average signal and calculate mean value
data_aux=mean(data,3);
mean_value = mean(mean(data_aux));
%Those values over 2*mean are considered brain matter and those below are
%set to 0 in the mask
for i=1:size(data,1)
    for j=1:size(data,1)
        if (mean_value*1.5) > data_aux(i,j)
            mask(i,j) = 0;
        else
            mask(i,j) = 1;
        end
    end
end

%The mask is eroded to avoid skull values still remaining
%mask = erodemask(mask,1);