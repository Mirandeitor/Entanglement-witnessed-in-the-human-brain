function [data,parameters] = read_philips_data(Options)
%
% Fucntion which is going to read the data and the header parameters from
% the Philips data scanner, depending on the sort of Philips data.
%
% [A,B] = read_philips_data(Options) Give back the Data A and the
% Parameters B from depending on the value of Options.protocol
% Protocol can have the values: PAR or Nifti
%
% V 1.0 Creation of the document by David Lopez 2.12.2011
% V 1.1 Bug Fix in the input parameters 14.12.2011 by David Lopez
% V 1.2 Multishot case also added 21.03.2012 by David Lopez
%
%
fprintf('Im reading Philips data \n');

if nargin==0 || isempty(Options.protocol)
    protocol = 'Nifti';
    filename = '';
    filenamePARData = '';
end

if nargin==1 && isempty(Options.protocol)
    protocol = 'Nifti';
else
    protocol = Options.protocol;
end

if nargin==1 && ~isfield(Options,'shot')
    shot = 'shot';
else
    shot = Options.shot;
end

if nargin==1 && ~isfield(Options,'filename')
    filename = '';
elseif strcmp(protocol,'Nifti')
    filename = Options.filename;
else
    filename = Options.filename;
    filenameRECData = Options.filenameREC;
end

switch protocol
    %Nifti case
    case 'Nifti'
        OptionsNifti.filename = filename;
        OptionsNifti.protocol = 'Nifti';
        [data, parameters]=read_nifti(OptionsNifti);
    %PAR/REC case    
    case 'PAR'
        OptionsHeader.filename = filename;
        parameters = read_philips_header(OptionsHeader);
        OptionsPAR.filename = filenameRECData;
        OptionsPAR.protocol = 'PAR';
        OptionsPAR.shot = shot;
        data = get_PAR_data(parameters,OptionsPAR);
    otherwise
        errordlg('Only Nifti and PAR protocols are supported');
end
end