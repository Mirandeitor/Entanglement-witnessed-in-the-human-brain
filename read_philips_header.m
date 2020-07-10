function param = read_philips_header(Options)
%
%
% Function which is going to read the header parameters of the REC file
% from the path set in Options.path
%
% P = read_philips_header(Options)
%
% V 1.0 Creation of the document by David Lopez 1.12.2011 
% V 1.1 Bug Fix Reading data properly from PAR file by David Lopez 4.1.2012
% V 1.2 Bug Fix Validation of the extension set to .PAR or .par by David Lopez
%       4.1.2012
%
%

if nargin == 0 || (nargin == 1 && isempty(Options.filename))
    [fname,pname] = uigetfile('*.par','Select the REC Philips file ');
    filename = strcat(pname,fname);
    [path,name,ext]= fileparts(filename);
    while ~strcmp(ext,'.PAR')
        error('The selected file is not a REC Philips file ');
    end
    cd(pname);    
    if ~filename
        return
    end
end

if nargin == 1 && ~isempty(Options.filename)
    switch exist(Options.filename)
    % If it's a file return the file
    case 2
       filename = Options.filename;
    % If it's a directory search for the files inside
       [path,name,ext]= fileparts(filename);
       if ~strcmp(ext,'.par') && ~strcmp(ext,'.PAR')==2
            error('The selected file is not a REC Philips file ');
       end
    case 7
       [path,name,ext]=fileparts();
       if  strcmp(ext,'.PAR')==2 || strcmp(ext,'.par')==2
           filename=[direc,'par'];
       else
           disp('ERROR: could not find the par file in the specified directory...');
           disp('...please select now!"');
           cd(filename)
           [fname,pname] = uigetfile('*.*', 'Open par Philips file');
           filename = strcat(pname,fname);
           if ~filename
               return
           end
       end
    % If it's nothing selected before or don't exist
    otherwise
       disp('ERROR: file or directory not found...');
       disp('...please select now!');
       [fname,pname] = uigetfile('*.*', 'Open par Philips file');
       filename = strcat(pname,fname);
       if ~filename
            return
       end   
    end
end

%%Load the file and read all the parameters
fid = fopen(filename);
if fid<0
    error('The file cannot be opened');
else
    param.name = filename;
end

mode = -1;
nHC=0;
nIC=0;
nSC=0;
loop = 1;
while(true)
    str=fgetl(fid);
    if ~ischar(str)
        break;
    end 
    if isempty(str)
        continue;
    end
    if strfind(str,'= DATA DESCRIPTION FILE =')
        mode = 0;
    end
    if strfind(str,'= GENERAL INFORMATION =')
        mode = 1;
    end
    if strfind(str,'= PIXEL VALUES =')
        mode = 2;
    end
    if strfind(str,'= IMAGE INFORMATION DEFINITION =')
        mode = 3;
    end
    if strfind(str,'= IMAGE INFORMATION =')
        mode = 4;
    end
    if strfind(str,'= END OF DATA DESCRIPTION FILE =')
        mode = 5;
    end
    %
    switch (mode)
        case -1;
        %Case 0 Reads the paramrmation about the description of the file
        case 0
        %Case 1 Reads the general paramrmation about the file
        case 1
            [type data] = get_general_information(str);
            switch(type)
                case 'Patientname'
                    param.(type)= data;
                case 'Protocolname'
                    param.(type)= data;
                case 'Examinationname'
                    param.(type)= data;
                case 'ExaminationdateTime'
                    param.(type)= data;
                case 'Seriestype'
                    param.(type)= data;
                case 'Acquisitionnr'
                    param.(type)= data;
                case 'Reconstructionnr'
                    param.(type)= data;
                case 'Scandurationsec'
                    param.(type)= data;
                case 'Maxnumberofcardiacphases'
                    param.(type)= data;
                case 'Maxnumberofechoes'
                    param.(type)= data;
                case 'MaxnumberofslicesLocations'
                    param.(type)= data;
                case 'Maxnumberofdynamics'
                    param.(type)= data;
                case 'Maxnumberofmixes'
                    param.(type)= data;
                case 'Patientposition'
                    param.(type)= data;
                case 'Preparationdirection'
                    param.(type)= data;
                case 'Technique'
                    param.(type)= data;
                case 'Scanresolutionxy'
                    param.(type)= data;
                case 'Scanmode'
                    param.(type)= data;
                case 'Repetitiontimems'
                    param.(type)= data;
                case 'Fovapfhrlmm'
                    param.(type)= data;
                case 'Waterfatshiftpixels'
                    param.(type)= data;
                case 'Angulationmidsliceapfhrldegr'
                    param.(type)= data;
                case 'Offcentremidsliceapfhrlmm'
                    param.(type)= data;
                case 'Flowcompensationnoyes'
                    param.(type)= data;
                case 'Presaturationnoyes'
                    param.(type)= data;
                case 'PhaseencodingvelocitycmSec'
                    param.(type)= data;
                case 'Mtcnoyes'
                    param.(type)= data;
                case 'Spirnoyes'
                    param.(type)= data;
                case 'Epifactornoepi'
                    param.(type)= data;
                case 'Dynamicscannoyes'
                    param.(type)= data;
                case 'Diffusionnoyes'
                    param.(type)= data;
                case 'Diffusionechotimems'
                    param.(type)= data;
                case 'Maxnumberofdiffusionvalues'
                    param.(type)= data;
                case 'Maxnumberofgradientorients'
                    param.(type)= data;
                case 'NumberoflabeltypesnoASL'
                    param.(type)= data;
                otherwise
                    %param.(type)= data;
            end
        %Case 2 Read paramrmation about the Pixel values (normally
        %paramrmation is in the REC file)
        case 2
        %Case 3 Image paramrmation definition    
        case 3
            [type datatype datalength] = get_image_paramrmation(str);
            if ~isempty(type)
                nIC=nIC+1;
                ImageparamrmationTags(nIC).Name=type;
                ImageparamrmationTags(nIC).DataType=datatype;
                ImageparamrmationTags(nIC).NumberOfValues=datalength;
            end
        %Case 4 Image parammration 
        case 4
            if(str(1)~='#');
                nSC=nSC+1;
                %Matches the expression
                vals=regexp(str, '\s+','split');
                vald=sscanf(str, '%lf')';
                current_loc=0;
                for i=1:length(ImageparamrmationTags)
                    IIT=ImageparamrmationTags(i);
                    if(strcmp(IIT.DataType,'string'))                      
                        Sliceparamrmation(nSC).(IIT.Name)=vals{current_loc+1};
                    else
                        Sliceparamrmation(nSC).(IIT.Name)=vald(current_loc+1:current_loc+IIT.NumberOfValues);
                    end
                    current_loc=current_loc+IIT.NumberOfValues;
                end
            end
        %Case 5 End of the file (Nothing is done)
        case 5
            break;
    end
    loop = loop + 1;
end
fclose(fid);

%param.HeaderComment=HeaderComment;

param.Sliceparamrmation=Sliceparamrmation;

param.ImageparamrmationTags=ImageparamrmationTags;

% Add Dimensions and Voxel Spacing. Warning, based on only 1 slice!
infof=param.Sliceparamrmation(1);
if(isfield(infof,'ReconResolution'))
    if(isfield(param,'MaxNumberOfSlicesLocations'))
        zs(1)=param.MaxNumberOfSlicesLocations;
        zs(2)=length(Sliceparamrmation)/zs(1);
        if((mod(zs(2),1)>0)||zs(2)==1)
            zs=length(Sliceparamrmation);
        end
    else
        zs=length(Sliceparamrmation);
    end
    
    param.Dimensions=[infof.ReconResolution zs];
else
    param.Dimensions=[param.ScanResolution length(Sliceparamrmation)];
end
if(isfield(infof,'PixelSpacing'))
    if(isfield(infof,'SliceThickness')&&isfield(infof,'SliceGap'))
        zs=infof.SliceThickness+infof.SliceGap;
    else
        zs=0;
    end
    param.Scales=[infof.PixelSpacing zs];
else
    param.Scales=[0 0 0];
end


[folder,filen]=fileparts(param.name);
param.FilenameREC=fullfile(folder,[filen '.rec']);
   
    
% Add bith depth
if(infof.ImagePixelSize)
    param.BitDepth=infof.ImagePixelSize;
else
    if(exist(param.FilenameREC,'file'))
        file_info=dir(param.FilenameREC);
        bytes=file_info.bytes;
        param.BitDepth=(bytes/prod(param.Dimensions))*8;
    end
end

end