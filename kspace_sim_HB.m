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
% - incorporate temporal data (scale one full beat by 'slow_factor')------done

% - detemine error interval during acquisition, ie what is the chance that
%   central kspace is filled. What is the 95% confidence interval for any
%   one acquisition (assuming numlines and slow factor are
%   constant).

% - make sure code is flexible to use experimentally determined
% values

% NOTES: 
% - the imaginary reconstruction looks to create a really nice angio image,
%   is this something worth investigating for angio images?
clear mrstruct_flow mrstruct_mag ;
bitdepth    = 2^12; % dicoms are encoded 12 bit
venc        = 60;   % cm/s (AJB: hardcoded this, obtained from header)

% SIM SETTINGS
numlines = 2;      % must be even, number of center lines to fill with slow data
slow_factor = 0.5; % change velocities in anomalous beat by this multiple
t = 5;  % time-frame
f = 1; % deviation from center f = (0,length of k-space], ie f = 1 is one line below center 

%[filename, filepath] = uigetfile({'*.dcm;*.ima'},'select all mag files','Multiselect','On');
%[filename1, filepath1] = uigetfile({'*.dcm;*.ima'},'select all flow files','Multiselect','On');
 filepath = uigetdir(pwd,'select mag file');
 filename = [dir(fullfile(filepath,'*.ima')); dir(fullfile(filepath,'*.dcm'))];
     filepath1 = uigetdir(pwd,'select flow file');
 filename1 = [dir(fullfile(filepath1,'*.ima')); dir(fullfile(filepath1,'*.dcm'))];

for i = 1:length(filename)
    y = filename(i).name;
    y1 = filename1(i).name;
[mrstruct_mag{i},~] = dicom_read_singlefile([filepath '\' y],1,[],0);
[mrstruct_flow{i},~] = dicom_read_singlefile([filepath1 '\' y1],1,[],0);
flow = mrstruct_flow{1,i}.dataAy;
mag = mrstruct_mag{1,i}.dataAy;
magstruct(:,:,i) = mag;
flowstruct(:,:,i) = flow;
end
bitrange = bitdepth./2;
mrstruct_flow = (flowstruct-bitrange)./bitrange; %result is +/- 1
mrstruct_flow = pi().*mrstruct_flow;



%% plot magnitude/phase along with reconstructed real/imaginary

% create imaginary and complex from magnitude and phase
mag   = magstruct;  % magnitude
phase = mrstruct_flow; % phase 
z = zeros(size(mag,1),size(mag,2),size(mag,3));
kspace_z = zeros(size(mag,1),size(mag,2),size(mag,3));
kspace_phase = zeros(size(mag,1),size(mag,2),size(mag,3));
for i = 1:size(mag,3)
z(:,:,i)= mag(:,:,i).*exp(1i*phase(:,:,i)); 
kspace_z(:,:,i)     = fft2(z(:,:,i));
kspace_phase(:,:,i) = fft2(phase(:,:,i));
end
% complex

% perform FFT on complex data (and phase to test if the same)
%kspace_z     = fft2(z);
%kspace_phase = fft2(phase);

% plot
hfig1   = figure('position',[1   405   1019   586]);
hax1 = subplot(3,3,1);
him1 = imagesc(mag(:,:,t));
title('magnitude')

hax2 = subplot(3,3,2);
him2 = imagesc(real(z(:,:,t)));
title('real(z)')

hax3  = subplot(3,3,3);
him3  = imagesc(real(ifft2(kspace_z(:,:,t))));
title('real(ifft2(kspace))')

hax4 = subplot(3,3,4);
him4 = imagesc(phase(:,:,t));
title('phase [-pi pi]')

hax5 = subplot(3,3,5);
him5 = imagesc(imag(z(:,:,t)));
title('imaginary(z)')

hax6  = subplot(3,3,6);
him6  = imagesc(imag(ifft2(kspace_z(:,:,t))));
title('imag(ifft2(kspace))')

hax7  = subplot(3,3,7);
him7  = imagesc(abs(fftshift(kspace_phase(:,:,t))));
caxis([0 1000])
title('fft(phase) | abs(kspace_{phase})')

hax8  = subplot(3,3,8);
him8  = imagesc(abs(fftshift(kspace_z(:,:,t))));
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
him1 = imagesc(phase(:,:,t));
title('phase_{norm} [-pi pi]')

% get roi and compute masks/phase_slow
roi_tube   = imellipse(gca,[88.6878751204264 69.4691651891066 15.4136815720394 16.1123707399367]);
mask       = createMask(roi_tube,him1);
masks = zeros(size(mag,1),size(mag,2),size(mag,3));
for i = 1:size(mag,3)
    masks(:,:,i) = mask;
end
mask_norm  = mask.*phase;             % masked phase (for 'flow' computation)
mask_slow  = mask_norm.*slow_factor; % masked slow phase (for 'flow' computation)


% seed with original data
for time = 1:size(mag,3)
    for j = 1:size(mag,2)
        for i = 1:size(mag,1)
            if mask(i,j) == 1
                phase_slow(i,j,time) = mask_slow(i,j,time);
            else
                phase_slow(i,j,time) = phase(i,j,time);
            end
        end
    end
end

%phase_slow(mask) = mask_slow(mask);

y = phase_slow(mask);

z_slow = zeros(size(mag,1),size(mag,2),size(mag,3));
kspace_z_slow = zeros(size(mag,1),size(mag,2),size(mag,3));
kspace_phase_slow = zeros(size(mag,1),size(mag,2),size(mag,3));
% compute complex slow data
for i = 1:size(mag,3)
z_slow(:,:,i) = mag(:,:,i).*exp(1i*phase_slow(:,:,i)); % complex data for slow phase/mag data

% perform FFT on complex data (and phase to test if this is the same)
kspace_z_slow(:,:,i)     = fft2(z_slow(:,:,i));
kspace_phase_slow(:,:,i) = fft2(phase_slow(:,:,i)); % won't play with this for now, doesn't seem worth it
end

hax2 = subplot(3,3,2);
him2 = imagesc(phase_slow(:,:,t));
title('phase_{slow} [-pi pi]')

hax3 = subplot(3,3,3);
him3 = imagesc(abs(fftshift(kspace_z(:,:,t))));
title('abs(kspace_z))')

hax4 = subplot(3,3,4);
him4 = imagesc(mask_norm(:,:,t));
title('mask_{norm}')
%compute mean phase shift across ROI
for x = 1:size(mag,3)
    
w = mask_norm(:,:,x); 
w(w == 0) = nan;
norm_mean(x) = nanmean(w(:)*60/pi());
end
w = mask_norm(:,:,t); 
w(w == 0) = nan;
vel_mean = nanmean(w(:))*60/pi();
text(81,90,['vel_{mean}=' num2str(vel_mean,3)])
hax5 = subplot(3,3,5);
him5 = imagesc(real(z_slow(:,:,t)));
title('real(z_{slow})')

hax6 = subplot(3,3,6);
him6 = imagesc(abs(fftshift(kspace_z_slow(:,:,t))));
title('abs(kspace_{zslow})')

hax7 = subplot(3,3,7);
him7 = imagesc(mask_slow(:,:,t));
title(['mask_{slow} ' num2str(slow_factor*100) '% norm'])
y = phase_slow.*mask;
y(y==0) = nan;
u = y(:,:,t);
%w = phase_slow(mask(:));
vel_slow_mean = nanmean(u(:))*60/pi();
text(81,90,['vel_{mean}=' num2str(vel_slow_mean,3)])

hax8 = subplot(3,3,8);
him8 = imagesc(imag(z_slow(:,:,t)));
title('real(z_{slow})')

hax9 = subplot(3,3,9);
him9 = imagesc(([abs(fftshift(kspace_z(:,:,t))) - abs(fftshift(kspace_z_slow(:,:,t)))]));
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
%mask_kcenter(1:numlines/2,:) = 1;
%mask_kcenter(103:104,:,:) = 1;
if round((size(mask_kcenter,1)-numlines/2)*(f/size(mask_kcenter,1))) == round(size(mask_kcenter,1)*(f/size(mask_kcenter,1)))
 mask_kcenter(round((size(mask_kcenter,1)-numlines/2)*(f/size(mask_kcenter,1))) :round(size(mask_kcenter,1)*(f/size(mask_kcenter,1)))+1,:,:) = 1;
else
mask_kcenter(round((size(mask_kcenter,1)-numlines/2)*(f/size(mask_kcenter,1))) :round(size(mask_kcenter,1)*(f/size(mask_kcenter,1))),:,:) = 1;
end
%create synthetic kspace with center containing slow data
kspace_z_syn = kspace_z;
for time = 1:size(mag,3)
    for i = 1:size(mag,1)
        for j = 1:size(mag,2)
            if mask_kcenter(i,j,time) == 1
                kspace_z_syn(i,j,time) = kspace_z_slow(i,j,time);
            else
                kspace_z_syn(i,j,time) = kspace_z(i,j,time);
            end
        end
    end
end
%kspace_z_syn(mask_kcenter(:)) = kspace_z_slow(mask_kcenter(:));

%convert to image data
z_syn     = ifft2(kspace_z_syn);
mag_syn   = abs(z_syn);
phase_syn = angle(z_syn);

syn_phase = kspace_z.*mask_kcenter;
% plot process of filling kspace
hfig3 = figure('position',[40   360   1019   586]);
hax1 = subplot(3,3,1);
him1 = imagesc(circshift(mask_kcenter(:,:,t),size(kspace_z,1)/2));
title(['center kspace mask, n=' num2str(numlines)])
text(8,8,['[' num2str(size(kspace_z)) ']'],'color','w')

%plot synthetic kspace with center containing slow data
hax2 = subplot(3,3,2);
him2 = imagesc(abs(fftshift(kspace_z_syn(:,:,t))));
title('synthetic kspace (center=slow)')

% show mag content of new synthesized data
hax3 = subplot(3,3,3);
him3 = imagesc(mag_syn(:,:,t));
title('synthetic image, mag')

% show image content of kspace center
hax4 = subplot(3,3,4);
him4 = imagesc(real(ifft2(kspace_z(:,:,t).*mask_kcenter(:,:,t))));
title('image content mask, real')

% show real image content of new synthesized data
hax5 = subplot(3,3,5);
him5 = imagesc(real(z_syn(:,:,t)));
title('synthetic image, real')

% show phase content of new synthesized data
hax6 = subplot(3,3,6);
him6 = imagesc(phase_syn(:,:,t));
title('synthetic image, phase')

% show image content of kspace center
hax7 = subplot(3,3,7);
him7 = imagesc(imag(ifft2(syn_phase(:,:,t))));
title('image content mask, imag')

% show imaginary image content of new synthesized data
hax8 = subplot(3,3,8);
him8 = imagesc(imag(z_syn(:,:,t)));
title('synthetic image, imag')

% show phase content of new synthesized data
hax9 = subplot(3,3,9);
him9 = imagesc(phase_syn(:,:,t));
title('synthetic image, phase')
for x = 1:size(mag,3)
    p = phase_syn(:,:,x).* mask;
    p(p==0) = nan;
    syn_mean(x) = nanmean(p(:))*60/pi();
end
p = phase_syn(:,:,t).*mask;
p(p==0) = nan;
vel_syn_mean = nanmean(p(:))*60/pi();
text(84,66,['vel_{mean}=' num2str(vel_syn_mean,3)])
% roi_tube   = imellipse(gca,[88.6878751204264 69.4691651891066 15.4136815720394 16.1123707399367]);

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
%%
q = 1;
while q < size(mag,1)+1
    mask_kcenter = false(size(kspace_z));
%mask_kcenter(1:numlines/2,:,:) = 1;
%f = .85;
if q == size(mag,1)
    mask_kcenter(1,:,:) = 1;
        mask_kcenter(q,:,:) = 1;
else
   mask_kcenter(q:q+1,:,:) = 1;
 
end


kspace_z_syn = kspace_z;
for time = 1:size(mag,3)
    for i = 1:size(mag,1)
        for j = 1:size(mag,2)
            if mask_kcenter(i,j,time) == 1
                kspace_z_syn(i,j,time) = kspace_z_slow(i,j,time);
            else
                kspace_z_syn(i,j,time) = kspace_z(i,j,time);
            end
        end
    end
end
z_syn     = ifft2(kspace_z_syn);
mag_syn   = abs(z_syn);
phase_syn = angle(z_syn);
for i = 1:size(mag,3)
    p = phase_syn(:,:,i).* mask;
    p(p==0) = nan;
    ty(q,i) = nanmean(p(:))*60/pi();
end

q = q+1;
end

%%
%for i = 1:size(mag,3)
[h,p,ci,~] = ttest(ty);
%con_int(:,i) = ci;
chance = (sum(ty < ci(1,:))/size(ty,1))*100;

error_across_time_and_numlines = abs(norm_mean - ty);
%%
figure;
plot(1:13,error_across_time_and_numlines(:,1:13));