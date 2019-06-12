function row=metadata(pathname,filename)
%% DESCRIPTION
% This fcn extracts out the metadata information from a Nikon image
% obtained from the Nikon AZ100 confocal microscope.
%
%% INPUT VARIABLES
% pathname: pathname of the image files
% 
% filename: filename cell variable of the image file(s)
%% OUTPUT VARIABLES
% row: cell array containing the organized and metadata information
% 
%% EXTRACT METADATA
row=cell(length(filename),12);

for dum=1:length(filename)
    I=import_tiff_stack([pathname,filename{dum}],1,'skip',1);
    try
        a1=I.info.ImageDescription;

        %laser (%)
        ii=strfind(a1,'Laser Power')+length('Laser Power}: ');
        b=a1(ii(1):ii(1)+5);
        ii2=strfind(b,newline);
        laser=str2double(b(1:ii2-1));% power (%)

        %gain (%)
        ii=strfind(a1,'Gain')+length('Gain}: ');
        b=a1(ii(1):ii(1)+4);
        ii2=strfind(b,' ');
        gain=str2double(b(1:ii2-1));% gain

        %excitation
        ii=strfind(a1,'ExcitationWavelength')+length('ExcitationWavelength="');
        b=a1(ii(1):ii(1)+4);
        ii2=strfind(b,'"');
        ex=str2double(b(1:ii2-1));% excitation wavelength in nm

        %emission
        ii=strfind(a1,'EmissionWavelength')+length('EmissionWavelength="');
        b=a1(ii(1):ii(1)+5);
        ii2=strfind(b,'"');
        em=str2double(b(1:ii2-1));% emission wavelength in nm

        %nominal magnification
        ii=strfind(a1,'NominalMagnification')+length('NominalMagnification="');
        b=a1(ii(1):ii(1)+4);
        ii2=strfind(b,'"');
        zoom=str2double(b(1:ii2-1));% zoom
    catch
        gain=[];
        laser=[];
        zoom=[];
        ex=[];
        em=[];
        disp('No Image Description field found!');
    end

    % nominal image resolution
    res=I.info.UnknownTags(2).Value;%res in um/px

    % x, y, and z position as reported by Nikon image aqcuisition program
    x=I.info.UnknownTags(5).Value(1);% x pos in um
    y=I.info.UnknownTags(5).Value(2);% y pos in um
    z=I.info.UnknownTags(5).Value(3);% z pos in um

    % Put everything into a cell array
    row(dum,:)=...
        {pathname,filename{dum},x,y,z,res,gain,laser,zoom,ex,em,true};
end