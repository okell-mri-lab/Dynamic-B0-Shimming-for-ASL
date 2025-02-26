% This function rescales and normalises an image based on a minimum and
% maximum percentile intensity.  RealOrAbs determines whether the real
% ('r') or absolute ('a') values are taken for the output image.
%
% function outIm = RescaleIm(inIm, MinIntensityPerc, MaxIntensityPerc, RealOrAbs, UseAbsoluteScaling)

function outIm = RescaleIm(Im, MinIntensityPerc, MaxIntensityPerc, RealOrAbs, UseAbsoluteScaling)
  

  if nargin < 2; MinIntensityPerc = 0; end
  if nargin < 3; MaxIntensityPerc = 100; end
  if nargin < 4; RealOrAbs = 'a'; end
  if nargin < 5; UseAbsoluteScaling = false; end
  
  % Take the real or absolute value
  if lower(RealOrAbs) == 'r'
    outIm = real(Im);
	  
    % Remove negative components
    outIm(outIm < 0) = 0;
	  
  else  % Use the absolute
    %outIm = abs(Im)/max(abs(Im(:)));
    outIm = abs(Im);
  end

  % Decompose into a magnitude image and a component vector
  Mag = sqrt(sum(outIm.^2,4)); 
  Vec = outIm./repmat(Mag,[1 1 1 size(outIm,4)]);
  
  % Correct for vectors where magnitude is zero
  Vec(repmat(Mag,[1 1 1 size(outIm,4)])==0)=0;
  
  % Find the minimum and maximum intensity in the magnitude image to
  % display
  if UseAbsoluteScaling
    MinIntensity = MinIntensityPerc;
    MaxIntensity = MaxIntensityPerc;
  else
    MinIntensity = prctile(Mag(:), MinIntensityPerc);
    MaxIntensity = prctile(Mag(:), MaxIntensityPerc);
  end
  
  % Rescale the magnitude image
  Mag = Mag - MinIntensity;
  MaxIntensity = MaxIntensity - MinIntensity;
  Mag(Mag < 0) = 0;
  Mag(Mag > MaxIntensity) = MaxIntensity;
  Mag = Mag / MaxIntensity;

  % Reintroduce the vector components
  outIm = Vec .* repmat(Mag,[1 1 1 size(outIm,4)]);

%  tmp = sqrt(sum(outIm.^2,4)); max(tmp(:))