function mask_roi=draw_mask(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%yang created, update:Jan 15 2022%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% img is two dimension, one slice
I_tmp1=abs(varargin{1});
I=mean(I_tmp1,3);

if (nargin < 2)
  thres = 0.4;
else
  thres  = varargin{2};
end

mask_roi=[];
% mask_ini=I>mean2(I)*thres;
mask_ini=I>max(I(:))*thres;

range_high=mean(I(mask_ini>0));

mask=false(size(mask_ini));
i=1;

f=figure;
set(gcf,'outerposition',get(0,'screensize'));
H = uicontrol(f,'Style', 'PushButton', ...
    'String', 'Break', ...
    'Callback', 'delete(gcbf)');
H.Position = [900 90 50 20];
while true
    
%     subplot(1,2,1),show_overlay(I,mask_ini,1,range_high);
%     subplot(1,2,2),
    show_overlay(I,mask,2,range_high);%update new ROI
    
    
    mask_roi_tmp1=roipoly();
    
    if (ishandle(H))
        mask_roi_tmp = mask_roi_tmp1;%.*mask_ini;
    else
        break;
    end
    
    if  any(mask_roi_tmp(:) > 0.5)           %determine if the ROI is null
        mask_roi(:,:,i)=mask_roi_tmp;
        mask=mask+mask_roi(:,:,i);
        %subplot(1,2,2), show_overlay(I,mask,2);
        fprintf('Have Drawn No.%d ROI\n', i);
        i=i+1;
    end
end
fprintf('Draw ROI ended\n');

function  show_overlay(img,mask,ind,range_high)
I = img;
mask=double(mask);

%choose color
if ind==1
    R=0.1;G=0.1;B=1;
elseif ind==2
    R=0;G=255;B=0;
end

color = cat(3, R*ones(size(I)), G*ones(size(I)), B*ones(size(I)));
mask_color=color.*repmat(mask,[1,1,3]);

imshow(I,[0 range_high]);
hold on
h = imshow(mask_color);
hold off

set(h, 'AlphaData', 0.3);
