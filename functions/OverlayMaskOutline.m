% Overlay the outline of a mask on the specified image using the given colour
%
% OutIm = OverlayMaskOutline(Im,Mask,col)

function OutIm = OverlayMaskOutline(Im,Mask,col,NormaliseImFirst,LineWidth,MaskLargerOrSmaller,ReplaceValues,UseMarkers,NewWin)

if nargin < 3 || isempty(col); col = [1 1 0]; end
if nargin < 4 || isempty(NormaliseImFirst); NormaliseImFirst = true; end
if nargin < 5 || isempty(LineWidth); LineWidth = 1; end
if nargin < 6 || isempty(MaskLargerOrSmaller); MaskLargerOrSmaller = 'larger'; end
if nargin < 7 || isempty(ReplaceValues); ReplaceValues = false; end
if nargin < 8 || isempty(UseMarkers); UseMarkers = false; end
if nargin < 9 || isempty(NewWin); NewWin = false; end

% Squeeze to ensure the third dimension is colour
if ndims(Im) >= 3 && (ndims(Im) == ndims(Mask))
    % Here we have a multi-slice/time image so concatenate
    Im = CatSlices(Im,true);
    Mask = CatSlices(Mask,true);
end
Im = squeeze(Im); Mask = squeeze(Mask);

% Normalise the image first if requested
if NormaliseImFirst
  Im = Im/max(Im(:));
end

% Check the third dimension has three colours
Sz = size(Im); if ndims(Im) < 3; Sz(3) = 1; end
if Sz(3) ~= 3
  if Sz(3) == 1  % Just pad to make a grayscale image
    Im = repmat(Im,[1 1 3]);
  else % Replace with VEPCASL colours
    Im = squeeze(Components2RGB(permute(Im,[1 2 4 3])));
  end
end

% Dilate/erode the mask
if strcmp(lower(MaskLargerOrSmaller),'larger')
  DilMask = imdilate(Mask,strel('disk',LineWidth)); 
else
  DilMask = imerode(Mask,strel('disk',LineWidth));   
end

% Remove the centre
OutlineMask = abs(DilMask - Mask); 

if NewWin; figure; end

if UseMarkers
    
    OutIm = Im;
    [x y] = meshgrid(1:size(Im,2),1:size(Im,1));
    xp = x(OutlineMask>0);
    yp = y(OutlineMask>0);
    
    imagesc(OutIm); axis equal; axis off;
    hold on;
    plot(xp,yp,'.','MarkerFaceColor',col,'MarkerSize',10);
    
else
    % Add the mask to the image
    if ReplaceValues
        for ii = 1:3
            tmp = Im(:,:,ii);
            tmp(OutlineMask>0) = col(ii);
            OutIm(:,:,ii) = tmp;
        end
    else
        % Set the colour
        FinalMask = zeros([size(OutlineMask) 3]);
        for ii = 1:3
            FinalMask(:,:,ii) = OutlineMask * col(ii);
        end
        
        OutIm = Im + FinalMask;
    end
    
    % Get rid of values above 1 and below 0
    OutIm(OutIm>1) = 1;
    OutIm(OutIm<0) = 0;

    % Display
    imagesc(OutIm); axis equal; axis off;
end


