function dataArray = dilatemask(dataArray,nDilate)


%DILATEMASK Dilates mask.
%
%   DILATEMASK Dilates a mask.
%
%   Syntax
%
%   DILATEMASK(A,n)
%
%
%
%   Description
%
%   DILATEMASK(A,n) dilates n times the N-D array A.
%
%
%   (Function supports in-place call using A.)
%
%
%   See also: ERODEMASK, IMDILATE

% F Schweser, 2009/07/23, mail@ferdinand-schweser.de




if nargin < 2
    error('Function requires two input arguments.')
end



dilateKernel = ones(3,3,3);
for jDilate = 1:nDilate
    dataArray = imdilate(dataArray,dilateKernel);
end



end