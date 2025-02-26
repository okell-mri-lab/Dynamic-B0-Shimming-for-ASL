% Calculate gradient offset and global frequency offset terms required to
% achieve good B0 homogeneity within the neck for PCASL labelling.
%
% [deltGx, deltGy, deltGz, Const_freqOff] = calc_pcasl_dyn_shim(rawFname, roi_exist, autoroi, r_mm, Frac, Shim2DOnly, GlobFreqCorrOnly)
%
% Inputs:
%   rawFname                Raw meas.dat file name of the dual-echo fieldmap
%                           acquisition in the neck
%   roi_exist               Flag to state that mask_shim.mat exists in the current
%                           directory/path and should be loaded rather than generating
%                           a new mask
%   autoroi                 Flag to state a new ROI should be defined automatically by
%                           the user selecting the centre of each vessel within the
%                           central slice
%   r_mm                    Radius of the initial cylindrical masks generated around each 
%                           vessel with autoroi
%   Frac                    Fraction of the 95th percentile magnitude signal that is used 
%                           to threshold the initial autoroi mask
%   Shim2DOnly              Apply a 2D correction only (using the central
%                           slice of the fieldmap if there are multiple slices).
%   GlobFreqCorrOnly        Apply a global frequency correction only (no
%                           linear terms)
%
% Outputs:
%   deltGx, deltGy, deltGz  Change in the linear shim terms required to
%                           correct for B0 inhomogeneity, in uT/m
%   Const_freqOff           Global frequency offset required for
%                           correction, in Hz

function [deltGx, deltGy, deltGz, Const_freqOff] = calc_pcasl_dyn_shim(rawFname, roi_exist, autoroi, r_mm, Frac, Shim2DOnly, GlobFreqCorrOnly)

warning('on','all')

% Find the directory of this file
mfiledir = fileparts(mfilename("fullpath"));

% Add the functions subdirectory
funcdir = fullfile(mfiledir,'functions');
disp(['Adding functions subdirectory to the path: ' funcdir])
addpath(funcdir);

% Set defaults for parameters not provided or empty
if nargin < 1 || isempty(rawFname);         rawFname = togetfile('Select raw meas.dat file');  end
if nargin < 2 || isempty(roi_exist);        roi_exist = false;                                 end
if nargin < 3 || isempty(autoroi);          autoroi = true;                                    end
if nargin < 4 || isempty(r_mm);             r_mm = 4;                                          end
if nargin < 5 || isempty(Frac);             Frac = 0.9;                                        end
if nargin < 6 || isempty(Shim2DOnly);       Shim2DOnly = false;                                end
if nargin < 7 || isempty(GlobFreqCorrOnly); GlobFreqCorrOnly = false;                          end

%% Load B0 data from twix
field_twix = mapVBVD(rawFname);

[img_M,img_Punw,dte,img_Pini] = get_fieldmap(field_twix);

if Shim2DOnly % Just select the central slice
    slcInd = ceil(size(img_M,3)/2);
else % Use all slices
    slcInd=1:size(img_M,3);
end

img_mag = img_M   (:,:,slcInd);
img_P1  = img_Punw(:,:,slcInd);
img_P2  = img_Pini(:,:,slcInd);

phz_ini = img_P1;%use unwrapped data
nx = size(img_mag,1);ny = size(img_mag,2);nz=size(img_mag,3);

para = readlocation(field_twix);
posX = para.positionX;%mm
posY = para.positionY;%mm
posZ = para.positionZ(slcInd);%mm

PX = repmat(reshape(posX,[nx 1 1]),[1 ny nz]);
PY = repmat(reshape(posY,[1 ny 1]),[nx 1 nz]);
PZ = repmat(reshape(posZ,[1 1 nz]),[nx ny 1]);
CS = ones(nx,ny,nz);
H  = [PX(:) PY(:) PZ(:) CS(:)];

%%  draw ROI
RotIndex=3;

if roi_exist % Read from file if requested
    load('mask_shim.mat')


elseif autoroi % More rapid ROI selection from single points
    VoxSz_mm = abs(posX(2)-posX(1)); % Assume isotropic
    r_vox = r_mm / VoxSz_mm;
    mask_m = rot90(point_mask(rot90(img_mag,RotIndex),r_vox,Frac),-RotIndex);
    save('mask_shim.mat',"mask_m")

else  % Standard ROI drawing for each slice separately    
    for i=1:nz
        mask_m0 = sum(draw_mask(rot90(img_mag(:,:,i),RotIndex),0.5),3);
        mask_m(:,:,i)=rot90(mask_m0,-RotIndex);

        figure;
        scale=200;
        subplot(121),
        imshow(rot90(mask_m(:,:,i),RotIndex),[]);colorbar;colormap 'jet';title('initial');
        subplot(122),
        aa=(img_Pini(:,:,i)-img_Punw(:,:,i)).*mask_m(:,:,i)/(2*pi*dte/1000);
        imshow(rot90(aa,RotIndex),[-scale scale]);colorbar;colormap 'jet';title('wrapped pixel');

        save('mask_shim.mat',"mask_m")
    end
end


%% calculate linear terms of the shimming gradient
phz_compen = -phz_ini(:);
mask_v = mask_m(:);
HWH = H'*((mask_v*ones(1,4)).*H);
aw  = pinv(HWH)*H'*(mask_v.*phz_compen);

rotMat=[1 0 0;0 1 0;0 0 1]; % AP direction
gamma = 42.58; % MHz/T
dGx =  aw(1)*1e6/(dte*gamma*2*pi);%μT/m
dGy =  aw(2)*1e6/(dte*gamma*2*pi);%μT/m
dGz =  aw(3)*1e6/(dte*gamma*2*pi);%μT/m

deltG  = rotMat*[dGx; dGy; dGz];
deltGx = deltG(1);%μT/m
deltGy = deltG(2);%μT/m
deltGz = deltG(3);%μT/m
Const_freqOff = aw(4)/(2*pi*dte/1000);% Hz

% For a single slice fieldmap, there is ambiguity between Gz and the global
% frequency offset, so force Gz to zero here and correct the global
% frequency term
if nz==1
    Const_freqOff = -(gamma*deltGz*posZ/1000+Const_freqOff);
    deltGz = 0;
else % Just correct the sign
    Const_freqOff = -Const_freqOff;
end

% global frequency offset of initial map
a1=phz_ini/(2*pi*dte/1000);
gf1=mean(a1(mask_m>0));

% For the global offset correction, set the linear terms to zero and update
% the frequency correction required
if GlobFreqCorrOnly
    deltGx = 0; deltGy = 0; deltGz = 0;
    Const_freqOff = gf1;
end

% Print the correction factors to the screen
if GlobFreqCorrOnly % Full 2D/3D dynamic shimming
    fprintf('\n%%%%%%%%%% Global frequency offset correction method %%%%%%%%%%\n')
elseif Shim2DOnly
    fprintf('\n%%%%%%%%%%  2D dynamic shimming method %%%%%%%%%%\n')
else
    fprintf('\n%%%%%%%%%%  3D dynamic shimming method %%%%%%%%%%\n')
end
fprintf('   Input the following values into the special card:\n')

fprintf('   X-shim = %.1f μT/m \n', deltGx);
fprintf('   Y-shim = %.1f μT/m \n', deltGy);
fprintf('   Z-shim = %.1f μT/m \n', deltGz);
fprintf('   FreqZ offset = %.1f Hz\n', Const_freqOff);

%% Simulated B0 maps
dGx_s= deltGx;  %uT/m
dGy_s= deltGy;  %uT/m
dGz_s= deltGz;  %uT/m
FreqOff_s=-Const_freqOff;% Hz % T.O. tweak to be consistent with above

aw_s(1,1)=dGx_s*(dte*gamma*2*pi)/1e6;
aw_s(2,1)=dGy_s*(dte*gamma*2*pi)/1e6;
aw_s(3,1)=dGz_s*(dte*gamma*2*pi)/1e6;
aw_s(4,1)=FreqOff_s*(2*pi*dte/1000);

phz_compen_s=reshape(H*aw_s,nx,ny,nz);
phz_sim=phz_compen_s+phz_ini;


%% Plot pre- and post-dynamic shim fieldmaps
scale=350;
ratio=2*pi*dte/1000;
for i=1:nz

    sim_p_m=phz_sim(:,:,i).*mask_m(:,:,i)/ratio;
    ini_p_m=phz_ini(:,:,i).*mask_m(:,:,i)/ratio;

    sim_p=phz_sim(:,:,i).*1/ratio;
    ini_p=phz_ini(:,:,i).*1/ratio;
    mask_rm=ini_p~=0;
    figure;
    warning('off','all')
    subplot(221),imshow(rot90(ini_p,3),[-scale scale]);colorbar;colormap 'jet'; title('before DynShim');
    subplot(222),imshow(rot90(abs(ini_p_m),3),[0 60]);colorbar;colormap 'jet'; title(['abs-initial   ',num2str(round(gf1)),' Hz']); 
    subplot(223),imshow(rot90(sim_p,3),[-scale scale]);colorbar;colormap 'jet'; title('after DynShim')
    subplot(224),imshow(rot90(abs(sim_p_m),3),[0 60]);colorbar;colormap 'jet'; title('abs-simulated')
end

% T.O. Add histogram of off-resonance before and after correction
figure; histogram(phz_ini(mask_m>0)/ratio); hold on
histogram(phz_sim(mask_m>0)/ratio); 
title 'Off-resonance frequencies within the mask'
xlabel 'Frequency/Hz'; ylabel 'Count'
legend({'Before dyn shim','After dyn shim'})

