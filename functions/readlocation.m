function para = readlocation( varargin )
twix_obj=varargin{1};
twix_obj = twix_obj(length(twix_obj));

sSliceArray=twix_obj.hdr.MeasYaps.sSliceArray.asSlice;

% When the value for these coordinate is very small, e.g.
% sPosition_dSag = 1.343452e-9 then the read_meas_dat function that readin the data would not
% recognize it and will leave the array empty so fix it here

for count = 1:length(sSliceArray(1:end))
    
    
    if ~isfield(sSliceArray{count},'sPosition')
        sSliceArray{count}.sPosition.dCor = 0;
        sSliceArray{count}.sPosition.dSag = 0;
        sSliceArray{count}.sPosition.dTra = 0;
    else
        if ~isfield(sSliceArray{count}.sPosition,'dCor')
            sSliceArray{count}.sPosition.dCor = 0;
        end
        
        if ~isfield(sSliceArray{count}.sPosition,'dSag')
            sSliceArray{count}.sPosition.dSag = 0;
        end
        
        if ~isfield(sSliceArray{count}.sPosition,'dTra')
            sSliceArray{count}.sPosition.dTra = 0;
        end
    end
end




%

for count = 1:length(sSliceArray(1:end))
    if ~isfield(sSliceArray{count}.sNormal,'dCor')
        sSliceArray{count}.sNormal.dCor = 0;
    end
    
    if ~isfield(sSliceArray{count}.sNormal,'dSag')
        sSliceArray{count}.sNormal.dSag = 0;
    end
    
    if ~isfield(sSliceArray{count}.sNormal,'dTra')
        sSliceArray{count}.sNormal.dTra = 0;
    end
end


NormalVec = [sSliceArray{1}.sNormal.dSag, sSliceArray{1}.sNormal.dCor, sSliceArray{1}.sNormal.dTra].';

for i=1:length(sSliceArray(1:end))
    Pos(i,1) = [sSliceArray{i}.sPosition.dSag].';
    Pos(i,2) = -[sSliceArray{i}.sPosition.dCor].';%correct it
    Pos(i,3) = -[sSliceArray{i}.sPosition.dTra].';%correct it
end

% 
para.SlicePos = Pos*NormalVec;%note the 'SlicePos' order is in ascending order. No need reorder,
                              % but need check if the data is in a consistent order. 
para.Pos = Pos;

para.lBaseResolution     =twix_obj.hdr.Phoenix.sKSpace.lBaseResolution;
para.lPhaseEncodingLines =twix_obj.hdr.Phoenix.sKSpace.lPhaseEncodingLines;
para.lPartitions         =twix_obj.hdr.Phoenix.sKSpace.lPartitions;
para.alDwellTime         =twix_obj.hdr.Phoenix.sRXSPEC.alDwellTime{1};%para.alDwellTime=twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1};
para.dPhaseFOV           =twix_obj.hdr.Phoenix.sSliceArray.asSlice{1, 1}.dPhaseFOV;
para.dReadoutFOV         =twix_obj.hdr.Phoenix.sSliceArray.asSlice{1, 1}.dReadoutFOV;
para.dThickness_3d       =twix_obj.hdr.Phoenix.sSliceArray.asSlice{1, 1}.dThickness;
para.dThickness_2d       =twix_obj.hdr.Phoenix.sSliceArray.asSlice{1, 1}.dThickness;
para.NColMeas            =twix_obj.hdr.Meas.NColMeas;
para.lImagesPerSlab      =twix_obj.hdr.Phoenix.sKSpace.lImagesPerSlab;

para.dPhaseRes           =para.dPhaseFOV   / para.lPhaseEncodingLines;
para.dFreqRes            =para.dReadoutFOV / para.lBaseResolution;

para.ADC_duration        =para.alDwellTime * para.NColMeas/1000;
para.ReadoutOSFactor = twix_obj.hdr.Dicom.flReadoutOSFactor;
para.BWperPixel          =round(1000000/(para.alDwellTime/1000*para.ReadoutOSFactor)/para.lBaseResolution);


%%%---------------- positionY for PSF   
disp(['X center: ', num2str(para.Pos(1,1)), ' mm   '])
disp(['Y center: ', num2str(para.Pos(1,2)), ' mm   '])
disp(['Z center: ', num2str(para.Pos(  round((size(sSliceArray,2)+1)/2) ,3)), ' mm   '])

Pos = para.Pos;% correct it here
for i=1:size(Pos,1)
    Position = Pos(i,:);     % fov position
    if nargin < 2%AP
        y_offset = para.dPhaseRes/2 + Position(2) ;
        x_offset = para.dFreqRes/2 + Position(1) ;
    else %RL
        y_offset = para.dPhaseRes/2 + Position(1) ;
        x_offset = para.dFreqRes/2 + Position(2) ;
    end
    
    centerY = (1 + para.lPhaseEncodingLines)/2;
    SlicePosY = para.dPhaseRes * (centerY - [1:para.lPhaseEncodingLines]) + y_offset;
    para.positionY=SlicePosY;%%%%note here unit is mm
    
    centerX = (1 + para.lBaseResolution)/2;
    SlicePosX = para.dFreqRes * (centerX - [1:para.lBaseResolution]) + x_offset;
    para.positionX=SlicePosX;%%%%note here unit is mm    
    
    %%%---------------- positionZ for PSF
    para.positionZ(i)=Position(3);%%%%note here unit is cm
    %%%%%%%%%%%%%%%%%%%%
end