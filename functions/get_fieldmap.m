function [img_mag,img_phz,dte,img_phz2] = get_fieldmap( twix )
%% Load & Reconstruct B0 data from twix
twix = twix(length(twix));
twix.image.flagRemoveOS = true;

kdata = permute(twix.image(),[1 3 2 8 5 4 6 7]); 
im0 = squeeze(fft2call(kdata));


%twix_acs1.hdr.Config.relSliceNumber
slcN=size(im0,5);
slc_order= twix.hdr.Meas.chronSliceIndices(1:slcN);

[~,reorder_index]=sort(slc_order);
im=im0(:,:,:,:,reorder_index);

im_complex = sum(im(:,:,:,2,:).*conj(im(:,:,:,1,:)),3);
im_complex = squeeze(im_complex);                 
img_mag=squeeze(mean(sos(im,3),4));
img_phz=sunwrap(im_complex);%unwrap
img_phz2=angle(im_complex);


TE1 = twix.hdr.Phoenix.alTE{1} / 1e3; TE2 = twix.hdr.Phoenix.alTE{2} / 1e3;
dte=TE2-TE1;%ms

% figure;
% RotIndex=3;
% ratio=2*pi*dte/1000;
% scale=490;
% subplot(131),imshow(rot90(img_phz/ratio,RotIndex),[-scale scale]);colorbar;colormap 'jet'
% title('unwrap (Hz)')
% subplot(132),imshow(rot90(img_phz2/ratio,RotIndex),[-scale scale]);colorbar;colormap 'jet'
% title('wrap (Hz)')
% subplot(133),imshow(rot90((img_phz2-img_phz)/ratio,RotIndex),[-scale scale]);colorbar;colormap 'jet'
% title('difference (Hz)')
