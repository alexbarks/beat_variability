% Purpose: this script simulates the k-space of magnitude and phase
% values obtained from 2D through plane MRI dicoms. The motivation is to
% determine the impact of an anomalous beat during an ECG-gated PCMRI
% acquisition. For proof of concept, the first effort looks at a pulsatile
% phantom during peak flow and simulates an anomalous beat by scaling the
% phase data by a 'slow_factor' and uses the synthetic kspace lines from
% this slow data to fill the central kspace of the original data. The
% number of central lines filled with this slow data is determined by
% 'numlines'.

% Alex Barker, Northwestern University 20171105

% TODO:
% - incorporate temporal data (scale one full beat by 'slow_factor')
% - detemine error interval during acquisition, ie what is the chance that
%   central kspace is filled. What is the 95% confidence interval for any
%   one acquisition (assuming numlines and slow factor are constant).
% - make sure code is flexible to use experimentally determined values

% NOTES: 
% - the imaginary reconstruction looks to create a really nice angio image,
%   is this something worth investigating for angio images?

%% read data & set up variables


% DICOM INFO
bitdepth    = 2^12; % dicoms are encoded 12 bit
venc        = 60;   % cm/s (AJB: hardcoded this, obtained from header)
i=0;
for numlines=0:2:40;
    for slow_factor=0.5:0.05:1;
i=i+1;
% SIM SETTINGS
% numlines = 2;      % must be even, number of center lines to fill with slow data
% slow_factor = 0.5; % change velocities in anomalous beat by this multiple

% Complete time resolved pulsatile phantom data is located on the server at:
% 10.254.136.37\data_imaging\cv_mri\misc\in_vitro_Phantoms_Emory\pulsatile\FLOW_VENC-SCOUT_60-120_PULSATILE_0025\
% 10.254.136.37\data_imaging\cv_mir\misc\in_vitro_Phantoms_Emory\pulsatile\FLOW_VENC-SCOUT_60-120_PULSATILE_P_0026\

% local paths to phantom data
root_mag  = 'L:\cv_mri\beat_variability\data\20150325_phantom_pulsatile\FLOW_VENC-SCOUT_60-120_PULSATILE_0025\';
root_flow = 'L:\cv_mri\beat_variability\data\20150325_phantom_pulsatile\FLOW_VENC-SCOUT_60-120_PULSATILE_P_0026\';

% specific files for 'peak flow' (by eye), will work on one time point for prototyping
file_mag  = 'IN-VITRO_PHANTOM_PULSATILE.MR.CARDIAC_NEW_CARDIAC_NEW11-6-15.0025.0006.2015.03.24.11.29.12.973288.43664752.IMA';
file_flow = 'IN-VITRO_PHANTOM_PULSATILE.MR.CARDIAC_NEW_CARDIAC_NEW11-6-15.0026.0006.2015.03.24.11.29.12.973288.43665007.IMA';

path_mag  = [root_mag file_mag];
path_flow = [root_flow file_flow];

[mrstruct_mag ,info_mag]  = dicom_read_singlefile(path_mag ,1,[],0);
[mrstruct_flow,info_flow] = dicom_read_singlefile(path_flow,1,[],0); % need to shift and convert to pi
									  
bitrange = bitdepth./2;
mrstruct_flow.dataAy = (mrstruct_flow.dataAy-bitrange)./bitrange; %result is +/- 1
mrstruct_flow.dataAy = pi().*mrstruct_flow.dataAy;

%% plot magnitude/phase along with reconstructed real/imaginary

% create imaginary and complex from magnitude and phase
mag   = mrstruct_mag.dataAy;  % magnitude
phase = mrstruct_flow.dataAy; % phase 
z     = mag.*exp(1i*phase);   % complex

% perform FFT on complex data (and phase to test if the same)
kspace_z     = fft2(z);
kspace_phase = fft2(phase);

% plot
hfig1   = figure('position',[1   405   1019   586]);
hax1 = subplot(3,3,1);
him1 = imagesc(mag);
title('magnitude')

hax2 = subplot(3,3,2);
him2 = imagesc(real(z));
title('real(z)')

hax3  = subplot(3,3,3);
him3  = imagesc(real(ifft2(kspace_z)));
title('real(ifft2(kspace))')

hax4 = subplot(3,3,4);
him4 = imagesc(phase);
title('phase [-pi pi]')

hax5 = subplot(3,3,5);
him5 = imagesc(imag(z));
title('imaginary(z)')

hax6  = subplot(3,3,6);
him6  = imagesc(imag(ifft2(kspace_z)));
title('imag(ifft2(kspace))')

hax7  = subplot(3,3,7);
him7  = imagesc(abs(fftshift(kspace_phase)));
caxis([0 1000])
title('fft(phase) | abs(kspace_{phase})')

hax8  = subplot(3,3,8);
him8  = imagesc(abs(fftshift(kspace_z)));
caxis([0 1.5e05])
title('fft(z) | abs(kspace_z)')

% fft(phase) is not the same as fft(z) - obviously
hax9  = subplot(3,3,9);
% imagesc([abs(fftshift(kspace_z))-abs(fftshift(kspace_phase))])
% caxis([0 1.5e05])
% title('abs(kspace_z)-abs(kspace_{phase}')

axis([hax1 hax2 hax3 hax4 hax5 hax6 hax7 hax8 hax9],'image')
axis([hax1 hax2 hax3 hax4 hax5 hax6 hax7 hax8 hax9],'off')
colormap(hfig1,gray)

% clean up
clear temp_flow root_mag root_flow file_mag file_flow path_mag path_flow bit_depth bit_range
clear hax1 hax2 hax3 hax4 hax5 hax6 hax7 hax8 hax9 hfig1 him1 him2 him3 him4

%% create 'slow' flow phantom images
% make region in tube some% of velocity values and then create synthetic
% k-space, ie kspace_slow

% plot process of ROI and computation of 'slow flow'
hfig2 = figure('position',[20   380   1019   586]);
hax1 = subplot(3,3,1);
him1 = imagesc(phase);
title('phase_{norm} [-pi pi]')

% get roi and compute masks/phase_slow
roi_tube   = imellipse(gca,[88.6878751204264 69.4691651891066 15.4136815720394 16.1123707399367]);
mask       = createMask(roi_tube,him1);
mask_norm  = mask.*phase;             % masked phase (for 'flow' computation)
mask_slow  = mask_norm.*slow_factor;  % masked slow phase (for 'flow' computation)
phase_slow = phase;                   % seed with original data
phase_slow(mask) = mask_slow(mask);

% compute complex slow data
z_slow = mag.*exp(1i*phase_slow); % complex data for slow phase/mag data

% perform FFT on complex data (and phase to test if this is the same)
kspace_z_slow     = fft2(z_slow);
kspace_phase_slow = fft2(phase_slow); % won't play with this for now, doesn't seem worth it

hax2 = subplot(3,3,2);
him2 = imagesc(phase_slow);
title('phase_{slow} [-pi pi]')

hax3 = subplot(3,3,3);
him3 = imagesc(abs(fftshift(kspace_z)));
title('abs(kspace_z))')

hax4 = subplot(3,3,4);
him4 = imagesc(mask_norm);
title('mask_{norm}')
%compute mean phase shift across ROI
vel_mean = mean(phase(mask(:)))*60/pi();
text(81,90,['vel_{mean}=' num2str(vel_mean,3)])

hax5 = subplot(3,3,5);
him5 = imagesc(real(z_slow));
title('real(z_{slow})')

hax6 = subplot(3,3,6);
him6 = imagesc(abs(fftshift(kspace_z_slow)));
title('abs(kspace_{zslow})')

hax7 = subplot(3,3,7);
him7 = imagesc(mask_slow);
title(['mask_{slow} ' num2str(slow_factor*100) '% norm'])
vel_slow_mean = mean(phase_slow(mask(:)))*60/pi();
text(81,90,['vel_{mean}=' num2str(vel_slow_mean,3)])

hax8 = subplot(3,3,8);
him8 = imagesc(imag(z_slow));
title('real(z_{slow})')

hax9 = subplot(3,3,9);
him9 = imagesc([abs(fftshift(kspace_z)) - abs(fftshift(kspace_z_slow))]);
title('subtraction of kspace''s')

% make all the axes the same
hobj = findobj(gcf,'type','axes');
axis(hobj,'image')
axis(hobj,'off')
colormap(hfig2,gray)

%set caxis for phase images the same
set([hax1,hax2,hax4,hax7],'clim',[-pi() pi()],'xlim',[75 115],'ylim',[60 100])
set([hax3,hax6,hax9],     'clim',[0 1.5e5])

clear hax1 hax2 hax3 hax4 hax5 hax6 hax7 hax8 hax9 hfig2 him1 him2 him3 him4 him5 him6 him7 him8 him9 hobj

%% now synthetically fill center of kspace_z with kspace_z_slow

% make mask
mask_kcenter = false(size(kspace_z));
mask_kcenter(1:numlines/2,:) = 1;
mask_kcenter((end-numlines/2):end,:) = 1;

%create synthetic kspace with center containing slow data
kspace_z_syn = kspace_z;
kspace_z_syn(mask_kcenter(:)) = kspace_z_slow(mask_kcenter(:));

%convert to image data
z_syn     = ifft2(kspace_z_syn);
mag_syn   = abs(z_syn);
phase_syn = angle(z_syn);

% plot process of filling kspace
hfig3 = figure('position',[40   360   1019   586]);
hax1 = subplot(3,3,1);
him1 = imagesc(circshift(mask_kcenter,size(kspace_z,1)/2));
title(['center kspace mask, n=' num2str(numlines)])
text(8,8,['[' num2str(size(kspace_z)) ']'],'color','w')

%plot synthetic kspace with center containing slow data
hax2 = subplot(3,3,2);
him2 = imagesc(abs(fftshift(kspace_z_syn)));
title('synthetic kspace (center=slow)')

% show mag content of new synthesized data
hax3 = subplot(3,3,3);
him3 = imagesc(mag_syn);
title('synthetic image, mag')

% show image content of kspace center
hax4 = subplot(3,3,4);
him4 = imagesc(real(ifft2(kspace_z.*mask_kcenter)));
title('image content mask, real')

% show real image content of new synthesized data
hax5 = subplot(3,3,5);
him5 = imagesc(real(z_syn));
title('synthetic image, real')

% show phase content of new synthesized data
hax6 = subplot(3,3,6);
him6 = imagesc(phase_syn);
title('synthetic image, phase')

% show image content of kspace center
hax7 = subplot(3,3,7);
him7 = imagesc(imag(ifft2(kspace_z.*mask_kcenter)));
title('image content mask, imag')

% show imaginary image content of new synthesized data
hax8 = subplot(3,3,8);
him8 = imagesc(imag(z_syn));
title('synthetic image, imag')

% show phase content of new synthesized data
hax9 = subplot(3,3,9);
him9 = imagesc(phase_syn);
title('synthetic image, phase')
vel_syn_mean = mean(phase_syn(mask(:)))*60/pi();
text(84,66,['vel_{mean}=' num2str(vel_syn_mean,3)])
roi_tube   = imellipse(gca,[88.6878751204264 69.4691651891066 15.4136815720394 16.1123707399367]);

% make all the axes the same
hobj = findobj(gcf,'type','axes');
axis(hobj,'image')
axis(hobj,'off')
colormap(hfig3,gray)

%set caxis for kspace images the same
set(hax2,'clim',[0 1.5e5])
%set caxis for phase images the same
set([hax6,hax9],'clim',[-pi() pi()])
set(hax9,'xlim',[75 115],'ylim',[60 100])
close all
c(i,1)=numlines;
c(i,2)=slow_factor;
c(i,3)=vel_syn_mean;

    end
end

