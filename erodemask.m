function dataArrayObj = erodemask(dataArrayObj,nErode,boolUseFrame,type)

%ERODEMASK Erodes mask.
%
%   ERODEMASK Erodes a mask.
%
%   Syntax
%
%   ERODEMASK(A,n)
%   ERODEMASK(A,n,true)
%   ERODEMASK(A,K)
%
%
%
%   Description
%
%   ERODEMASK(A,n) erodes n times the N-D array A.
%
%   ERODEMASK(A,n,true) avoids artifacts at the edges (if trues in mask extend to
%   the edge of the matrix).
%
%   ERODEMASK(A,n,true, type) If type is equal to '2d' performs a 2d erosion if not it will do it in 3D

%   ERODEMASK(A,K) erodes the N-D array A using the kernel in K. 
%
%   
%
%   (Function supports in-place call using A.)
%
%
%   See also: DILATEMASK, IMERODE

% F Schweser, 2009/07/23, mail@ferdinand-schweser.de
%
%   v1.1 - 2009/12/15 - F Schweser - anti-aliasing frame 
%   v1.2 - 2009/02/15 - F Schweser - Now returns logical.
%   v1.3 - 2010/10/04 - F Schweser - Now works with mids object.
%   v1.4 - 2011/02/22 - A Deistung - BugFix: mids-output
%   v1.5 - 2011/02/28 - F Schweser - Kernel functionality added
%   v1.6 - 2011/05/31 - A Deistung - data input is now also supported if 
%                                    input is a struct with img field
%   v1.7 - 2011/06/30 - D Lopez    - If type is equal to 2d it erodes the mask in x
%                                    and y coordinates slice by slice. 2D
%                                    erosion.
%   v1.7.1 2011/06/30 - D Lopez    - Description updated
%   v1.7.2 2011/07/12 - D Lopez    - Type default value set up to '3d'



if nargin < 2
    error('Function requires two input arguments.')
end


if nargin > 2 && ~isempty(boolUseFrame)
    if ~islogical(boolUseFrame)
        error('Third parameter must be logical.')
    end
else
    boolUseFrame = false;
end

if nargin < 4
    %default mode set up to brain
    type = '3d';
end

if numel(nErode) == 1
    isKernel = false;
else
    isKernel = true;
    nErode = abs(nErode) / sum(abs(nErode(:))); % normalization
end
    




if isobject(dataArrayObj)
    dataArray = dataArrayObj.img;
elseif isstruct(dataArrayObj)
    if isfield(dataArrayObj, 'img')
        dataArray = dataArrayObj.img;
    else
        error('The img field is not available. Please pass a mids-object, struct with the field img, or a data array to erodemask!')
    end
else
    dataArray = dataArrayObj;
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% Add frame (if desired)
%%%%%%%%%%%%%%%%%%%%%%%%%
if boolUseFrame
    [dataArray unDovector] = prepareforconvolution(dataArray,[],2*(nErode+2));
end



if isKernel
    dataArray = ifftsave(fftsave(dataArray) .* fftsave(nErode),[],'symmetric');
    dataArray = dataArray > 0.95;
else
    if ~strcmp(type,'2d')
        erodeKernel = ones(3,3,3);
    else
        erodeKernel = ones(3,3);
    end
    for jErode = 1:nErode
        if ~strcmp(type,'2d')
            dataArray = imerode(dataArray,erodeKernel);
        else
            for jSlice = 1:size(dataArray,3)
                dataArray(:,:,jSlice) = imerode(dataArray(:,:,jSlice),erodeKernel);
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove frame
%%%%%%%%%%%%%%%%%%%%%%%%%
if boolUseFrame
    dataArray = prepareforconvolution(dataArray,[],[],unDovector);
end


dataArray = logical(dataArray);


if isobject(dataArrayObj) || isstruct(dataArrayObj)
    dataArrayObj.img = dataArray;
else 
    dataArrayObj = dataArray;
end


end