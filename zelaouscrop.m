function [ElementObj CropArray] = zealouscrop(ElementObj, cropArray, frameValue, bUndo)

%ZEALOUSCROP
%
%   Syntax
%
%   X = ZEALOUSCROP(A)
%   X = ZEALOUSCROP(A,c)
%   X = ZEALOUSCROP(A,[],v)
%   [X c] = ZEALOUSCROP(...)
%
%
%   Description
%
%   X = ZEALOUSCROP(A) removes all homogeneous slices, i.e. slices where all
%   elements are equal, from array A.
%
%   X = ZEALOUSCROP(A,c) crops A with information supplied by vector c. The
%   vector may be obtained from the function itself by using two output
%   variables.
%
%   X = ZEALOUSCROP(A,[],v) only slices where all values are equal to v
%   (real scalar) are removed. This makes the output of the function more
%   reliable in some extreme cases.
%
%   X = ZEALOUSCROP(A,c, [], bUndo) the cropping of the cropped matrix A
%   will be undone by the information supplied by vector c and setting
%   bUndo to true. The values in the added frames are set to zero.
%
%   X = ZEALOUSCROP(A,c, v, bUndo) the cropping of the cropped matrix A
%   will be undone by the information supplied by vector c and setting
%   bUndo to true. The values in the added frames are set to the value
%   given by v.

%   [X c] = ZEALOUSCROP(...) also returns a vector c that contains cropping
%   information and may be used in successive calls to obtain equal
%   results.


%% check input arguments
if nargin < 1 || isempty(ElementObj),      error('Not enough input arguments.');                       end


if isobject(ElementObj)
    Element = ElementObj.img;
else
    Element = ElementObj;
end



if islogical(Element)
    isInputLogical = true;
else
    isInputLogical = false;
end

if nargin < 4 || isempty(bUndo)
    bUndo = false;
end


if nargin < 2 || isempty(cropArray)
    
    if nargin < 3 || isempty(frameValue)
        frameValue = [];
    end
    
    %% init
    cropValue = [];
    
    dimensions = ndims(Element);
    
    first = zeros(dimensions, 1);
    last  = zeros(dimensions, 1);
    lastbeforecropping = zeros(dimensions, 1);
    sises = size(Element);
    lastbeforecropping(:) = sises;
    output = sprintf('return Value: %s(', inputname(1));
    
    
    %% main part
    for dim = 1 : dimensions
        for i = 1 : sises(1)
            %     if ~Test(Element(i, :), cropValue), break;  end
            if ~isUnique(Element(i, :), cropValue), break;  end
        end
        first(dim) = i;
        if dim == 1
            cropValue = unique(Element(i-1, :));
        end
        for i = sises(1) : -1 : first(dim)
            %if ~Test(Element(i, :)), break;  end
            if ~isUnique(Element(i, :), cropValue), break;  end
        end
        
        last(dim) = i;
        Element = Element(first(dim):last(dim), :);
        sises(1) = last(dim) - first(dim) + 1;
        Element = reshape(Element, sises);
        
        Element = shiftdim(Element, 1);
        sises = circshift(sises, [0 -1]);
        
        output = sprintf('%s%i:%i:%i, ', output, first(dim), last(dim), lastbeforecropping(dim));
    end
    
    if nargout > 1
        CropArray = [first, last, lastbeforecropping];
    else
        display(sprintf('%s)', output(1:numel(output)-2)));
    end
    
    
elseif ~bUndo
    % cropArray supplied and bUndo is false
    if ndims(Element) == 4
        ElementOut = zeros(cropArray(1,2)-cropArray(1,1)+1,cropArray(2,2)-cropArray(2,1)+1,cropArray(3,2)-cropArray(3,1)+1,size(Element,4));
        for jDim = 1:size(Element,4)
            ElementOut(:,:,:,jDim) = Element(cropArray(1,1):cropArray(1,2),cropArray(2,1):cropArray(2,2),cropArray(3,1):cropArray(3,2),jDim);
        end
        Element = ElementOut;
    elseif ndims(Element) == 3
        Element = Element(cropArray(1,1):cropArray(1,2),cropArray(2,1):cropArray(2,2),cropArray(3,1):cropArray(3,2));
    elseif ndims(Element) == 2
        Element = Element(cropArray(1,1):cropArray(1,2),cropArray(2,1):cropArray(2,2));
    else
        error('This number of dimensions is currently not supported.')
    end
    CropArray = cropArray;
    
else
    % cropArray supplied and bUndo is true
    if isempty(frameValue)
        frameValue = 0;
    end
    
    if ndims(Element) == 4
        ElementOut = zeros([cropArray(:,3)', size(Element, 4)]) + frameValue;
        ElementOut(cropArray(1,1):cropArray(1,2),cropArray(2,1):cropArray(2,2),cropArray(3,1):cropArray(3,2),:) = Element;
        
    elseif ndims(Element) == 3
        ElementOut = zeros(cropArray(:,3)')  + frameValue;
        ElementOut(cropArray(1,1):cropArray(1,2),cropArray(2,1):cropArray(2,2),cropArray(3,1):cropArray(3,2)) = Element;
    elseif ndims(Element) == 2
        ElementOut = zeros(cropArray(:,3)')  + frameValue;
        ElementOut(cropArray(1,1):cropArray(1,2),cropArray(2,1):cropArray(2,2)) = Element;
    else
        error('This number of dimensions is currently not supported.')
    end
    
    Element = ElementOut;
end



if isInputLogical
    Element = logical(Element);
end


if isobject(ElementObj)
    ElementObj.img = Element;
else
    ElementObj = Element;
end



%% member functions
    function isTrue = isUnique(Element, cropValue)
        uniqueVector = unique(Element);
        
        isTrue = numel(uniqueVector) == 1;
        
        if isTrue
            if uniqueVector ~= frameValue
                isTrue = false;
            end
            if uniqueVector ~= cropValue
                isTrue = false;
            end
        end
    end

end
