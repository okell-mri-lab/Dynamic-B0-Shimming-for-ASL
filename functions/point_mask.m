function mask_roi=point_mask(Im, r, Frac)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Tom Okell created, Feb 2025 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ImForDisp=mean(abs(Im),3); % Average over slices/echos if provided
ImForDisp=abs(Im(:,:,ceil(end/2))); % Use central slice

DispIm(ImForDisp,0,99); colormap gray;
title 'Click to select the centre of each vessel. Press "q" to finish.'
hold on;

% Initialize variables
points = []; % Store selected points
f = gcf;
set(f,'outerposition',get(0,'screensize'));

% Interactive point selection loop
while true
    [x, y, button] = ginput(1); % Get user click
 
    if button == 'q'
        disp('q pressed')
        break
    else
        points = [points; x, y]; % Store point
        plot(x, y, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Mark the point
    end
end

disp('Point selection finished.');
disp('Selected points:');
disp(points);

% Initialise output in 2D
mask_roi = false(size(ImForDisp));

% Define positions on a grid
[xg,yg] = meshgrid(1:size(Im,2),1:size(Im,1));

% Draw an ROI around each point
for ii = 1:size(points,1)
    mask_roi(((xg - points(ii,1)).^2 + (yg - points(ii,2)).^2) <= r ) = true;
end

% Extend across slices
mask_roi = repmat(mask_roi,[1 1 size(Im,3)]); 

% Show the intial mask
DispIm(mask_roi + abs(Im)/max(abs(Im(:))),0,99);
title 'Initial mask'

% Threshold at the relevant fraction of the 95th percentile
Thr = prctile(abs(Im(:)),95) * Frac;
mask_roi = mask_roi .* (abs(Im) >= Thr);

% Show the final mask
DispIm(mask_roi + abs(Im)/max(abs(Im(:))),0,99);
title 'Final mask'

end
