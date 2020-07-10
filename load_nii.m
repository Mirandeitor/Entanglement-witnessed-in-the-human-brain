function [data,param] = load_nii(parameters)

%
%Function which is going to return a N-dimensional array which the
%information about the image.
%
% A = load_nii(B) The N-Dimensional array A will be returned in function of
% the parameters of the header.
% [A,B]=load_nii() The N-Dimensional array A and the parameters B are
% returned
% *.img Not supported yet
%
% V 1.0 Creation of the document by David Lopez 13/12/2011
% V 1.1 Bug Fix by David Lopez 14/12/2011
% (V 1.X Read from .img files by David Lopez)
%

if nargin<1 || ~isstruct(parameters)
    param = read_nifti_header();
else
    param = parameters;
    
end

if ~isfield(param,'end')
    nii_end = true;
else
    nii_end = param.end;                
end


mode = 0;
counter = 0;
filename = param.Filename;
[path,name,ext] = fileparts(filename);
files =dir(path);
if strcmp(ext,'.nii')|| size(files,1) == 1
    mode = 1;
elseif strcmp(ext,'.nii')|| ~(size(files,1) == 1)
    [fname,pname] = uigetfile('*.*','Select the nii or img Philips file with the information of the image');
    filename = strcat(pname,fname);
    [path,name,ext] = fileparts(param.Filename);
    if strcmp(ext,'.nii')
        for j = 1:size(files,1)
            if (files(j).name == strcat(name,'.',ext))
                bytes = files(j).bytes; 
            end
        end
        for i = 1:size(files,1) 
            
            if (files(i).bytes == bytes)
                counter = counter +1;
            end
        end
        if counter > 1
            mode = 1;
        else
            mode = 2;   
        end
    elseif strcmp(ext,'.img')
        mode = 3;
    else 
        mode = 0;
    end
else
    error('There is only one file in this folder and it has not .nii extension');
end

fid = fopen(filename);
cd(path);
if fid < 0
    fprintf('could not open file %s\n',filename);
    return;
else
    switch mode
        % Mode = 1 --> everything is in the same file
        case 1
        	% Seek volume data start
            datasize=prod(param.Dimensions)*(param.BitVoxel/8);
            fsize=param.Filesize;
            % Search for the position of the information of the file
            fseek(fid,fsize-datasize,'bof');
            % Read Volume data
            switch(param.DataTypeStr)
                case 'INT8'
                    V = int8(fread(fid,datasize,'int8'));
                case 'UINT8'
                    V = uint8(fread(fid,datasize,'uint8'));
                case 'INT16'
                    V = int16(fread(fid,datasize,'int16'));
                case 'UINT16'
                    V = uint16(fread(fid,datasize,'uint16'));
                case 'INT32'
                    V = int32(fread(fid,datasize,'int32'));
                case 'UINT32'   
                    V = uint32(fread(fid,datasize,'uint32'));
                case 'INT64'
                    V = int64(fread(fid,datasize,'int64'));
                case 'UINT64'
                    V = uint64(fread(fid,datasize,'uint64'));    
                otherwise
                    V = uint8(fread(fid,datasize,'uint8'));
            end
            
        % Mode = 2 --> one or more separate files
        case 2    
        while(nii_end)
            counter = 0;
            datasize=prod(param.Dimensions)*(param.BitVoxel/8);
            files =dir(path);
            for i=1:size(files,1)
                filename = strcat(path,file(counter+i));
                [path, file, ext] = fileparts(filename);
                if (files(i).bytes == bytes)
                    if ~strcmp(ext,'nii')
                        switch(param.DataTypeStr)
                        case 'INT8'
                            V_par = int8(fread(fid,datasize,'int8'));
                        case 'UINT8'
                            V_par = uint8(fread(fid,datasize,'uint8'));
                        case 'INT16'
                            V_par = int16(fread(fid,datasize,'int16'));
                        case 'UINT16'
                            V_par = uint16(fread(fid,datasize,'uint16'));
                        case 'INT32'
                            V_par = int32(fread(fid,datasize,'int32'));
                        case 'UINT32'   
                            V_par = uint32(fread(fid,datasize,'uint32'));
                        case 'INT64'
                            V_par = int64(fread(fid,datasize,'int64'));
                        case 'UINT64'
                            V_par = uint64(fread(fid,datasize,'uint64'));      
                        otherwise
                            V_par = uint8(fread(fid,datasize,'uint8'));
                        end
                    end
                    V = V  + Vpar;
                end
            end
        end
        % This case is for the .img file 
        case 3
        otherwise
            fprint('It is not a supported data');
            return;
    end
end

fclose(fid);
% Reshape the volume data to the right dimensions
V = reshape(V,param.Dimensions);

data = V;
    
end