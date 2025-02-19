clc;clear;close all;warning('on','all')
if ismac || isunix
    addpath './functions'
else
    addpath '.\functions'
end

roi_exist=0;
autoroi = true; % T.O. Add option to define ROIs in a more automated way
r = 5; % mm, radius of initial masks

%% Load B0 data from twix
%field_twix = mapVBVD('meas_MID00098_FID09057_gre_field_mapping_xyz0.dat');
%field_twix = mapVBVD('~/Data/7T_ASL_Tests/2025-01-23/Raw_data/meas_MID00067_FID39515_to_gre_field_mapping_dyn_shim_online.dat');
%field_twix = mapVBVD('~/Data/7T_ASL_Tests/2025-02-04-phantom/Raw_data/meas_MID00101_FID40728_to_gre_field_mapping_dyn_shim_online.dat');
field_twix = mapVBVD('/Users/tokell/Data/7T_ASL_Tests/2025-02-19-subject/Raw_data/meas_MID00078_FID43014_yj_fieldmap_2mm_5slc.dat');
[img_M,img_Punw,dte,img_Pini] = get_fieldmap(field_twix);


slcInd=[1:size(img_M,3)];
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

if roi_exist
    load('mask_shim.mat')

    %     figure;
    %     ratio=2*pi*dte/1000;
    %     scale=350;
    %     imshow(rot90((img_phz2-img_phz1).*mask_m/ratio,RotIndex),[-scale scale]);colorbar;colormap 'jet'

elseif autoroi % More rapid ROI selection from single points
    mask_m = rot90(point_mask(rot90(img_mag,RotIndex),10,0.9),-RotIndex);
    save('mask_shim.mat',"mask_m")

else  % Standard ROI drawing
    %     mask_m0 = sum(draw_mask(rot90(img_tof(:,:,slcInd),RotIndex),0.5),3);
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


%% calculate linear term of the shimming gradient
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

fprintf('\n%%%%%%%%%%  2d&3d dynamic shimming method %%%%%%%%%%\n')
fprintf('   Input the following values into the special card:\n')
% global frequency offset of initial map
a1=phz_ini/(2*pi*dte/1000);
gf1=mean(a1(mask_m>0));

%fprintf('Global_FreqOff before DynShim_XY = %.1f Hz\n', gf1);
fprintf('   X-shim = %.1f μT/m \n', deltGx);
fprintf('   Y-shim = %.1f μT/m \n', deltGy);
if nz==1
    fprintf('   Z-shim = 0 μT/m \n');
    gf2 = -(gamma*deltGz*posZ/1000+Const_freqOff);
    fprintf('   FreqZ offset = %.1f Hz\n', gf2);   
else
    fprintf('   Z-shim = %.1f μT/m \n', deltGz);
    gf2 = -Const_freqOff;
    fprintf('   FreqZ offset = %.1f Hz\n', gf2);  
end
fprintf('\n%%%%%%%%%% Global frequency offset correction method %%%%%%%%%%\n')
fprintf('   Input the following values into the special card:\n')
fprintf('   X-shim = 0 μT/m \n');
fprintf('   Y-shim = 0 μT/m \n');
fprintf('   Z-shim = 0 μT/m \n');
fprintf('   FreqZ offset = %.1f Hz\n', gf1);

%% 
% simulated B0 maps
dGx_s= deltG(1);  %uT/m
dGy_s= deltG(2);  %uT/m
dGz_s= deltG(3); %uT/m
FreqOff_s=aw(4)/(2*pi*dte/1000);% Hz % T.O. tweak to be consistent with above

aw_s(1,1)=dGx_s*(dte*gamma*2*pi)/1e6;
aw_s(2,1)=dGy_s*(dte*gamma*2*pi)/1e6;
aw_s(3,1)=dGz_s*(dte*gamma*2*pi)/1e6;
aw_s(4,1)=FreqOff_s*(2*pi*dte/1000);

phz_compen_s=reshape(H*aw_s,nx,ny,nz);
phz_sim=phz_compen_s+phz_ini;

if nz<2
    aw_s(3,1)=0;
end
aw_s(4,1)=0;
phz_sim2=reshape(H*aw_s,nx,ny,nz)+phz_ini;
a2=(phz_sim2)/(2*pi*dte/1000);
m2=mean(a2(mask_m>0));
% fprintf('\n---Global_FreqOff after DynShim_XY = %.1f Hz\n', m2);

%% 
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

