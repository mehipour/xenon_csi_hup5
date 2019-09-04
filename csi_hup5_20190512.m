% Name: FA_Analysis_07_25_2016_FID30265_R80.m
% Purpose: Analyze spectroscopy data acquired with kr_C13_Calib R80
% Programmed: 07/26/2016
%
% update (1) by MP on 8/18/2016 for carbon CSI

% update (2) by MP on 5/4/2019 for xenon CSI

% w = warning('query','last');
% id = w.identifier;
% warning('off',id)

%% Data
% Control rat 1 Parameters

% Rat 60
FileNameRoot = '/Users/mehipour/Library/Mobile Documents/com~apple~CloudDocs/Data/HUP-5-new/20190510_Rat66/meas_MID00223_FID12783_kr_CSI_R143_20190501Cartesian_TE1_8';
%proton_path = '/Users/mehipour/Library/Mobile Documents/com~apple~CloudDocs/Data/HUP-5-new/20190510_Rat66/RAT_66.MR._.0052.0010.2019.05.10.18.14.06.486778.21827354.IMA';
ap_shift = 0;
fov = 300;
Nyshift = -2;
Nxshift = -6;

%% Input/Ouput Parameters
% visualization parameters
show_proton = 0;       % show proton image (gas phase)
zoom_y = 0.8;          % spectra y-scale
color_peaks = 0;     % make the integrated area of each peak colored.
spectrascale = 1;   % generally set to 1 to zoom in on the pyr region.
xgridshift = 0;      % shift mrsi along x
ygridshift = 0;      % shift mrsi along y

% acquisition parameters
Bandwidth = 1e9/122100;   % 1/(Dwell time)
Freq = 17612680;
centric = 0;
Nx = 24;
Ny = 24;
Np = 256;

% reconstruction parameters
show_real_mrsi = 0;      % if 1 then phase spectra
TE = -0.0002;            % echo time in seconds.     
baseline_poly_order = 0; % Baseline correction polynomial order 
lb = 100;                 % line broadening (Hz)

% Close all open figures
close all

%% Read raw data
InputFileName = [FileNameRoot, '.dat'];
twix = mapVBVD(InputFileName);

%% Rearrange Raw Data to make KSpace
% initialize variables 
% KSpaceArray = zeros(Nx,Ny,Np);
% raw_data = [];

% initialize additional matrices
n = 1:Np;       % sample vector
kspace = zeros(Nx,Ny,Np);
raw_data = [];

% Rearrange raw data into matrix
for ii=1:Nx
    %                  Col Cha Lin Par Sli Ave Phs Eco Rep Set Seg
    data = twix.image(  :,  1,  :,  1,  1,  1,  1,  1,  ii,  1,  1);
    raw_data = [raw_data data(:)];
    for jj = 1:Ny
        KSpaceArray(jj,ii,:) = squeeze(data(:,1,jj)).*exp(2i*pi*(ap_shift/fov*(jj-1))).*exp(-lb/Bandwidth*n');    
    end
end

% if proton (gas)background  image is availabe
if show_proton 
    proton = dicomread(proton_path);
else
    proton = zeros(Nx,Ny);
end

% rearrange spatial data in the k-space using the approrpiate trajectory
if centric
   CentricRO_Sampling_11_21_2016();
   gx = m_alTabXOrdering;
   gy = m_alTabYOrdering;
    for ii = 1:Ny
        for jj = 1:Nx
            kspace(gx(ii,jj),gy(ii,jj),:) = squeeze(KSpaceArray(jj,ii,:)); %.*exp(-n/ilb);
        end
    end
else
    kspace = KSpaceArray;
end

%% MRSI Reconstruction
% Show kspace
% figure(510);
% imagesc(abs(kspace(:,:,1)));

% Fourier transform
complex_img = rot90(fftshift(fftn((kspace))),1);
complex_img = circshift(complex_img,[Nyshift Nxshift 0]);
% if centric
complex_img = complex_img(:,:,Np:-1:1);
% end

% Find chemical shifts
find_xe129_chemical_shifts_hup5_20190505;

% Phase data
if ~show_real_mrsi
    mrsi = abs(complex_img);
    for ii = 1:Nx
        for jj = 1:Ny
           % second order baseine correction for magntitude spectra
           baseline = polyval(polyfit(baseline_idx,squeeze(mrsi(ii,jj,baseline_idx))',baseline_poly_order),n);
           mrsi(ii,jj,:) = squeeze(mrsi(ii,jj,:))-baseline';
        end
    end
else
    mrsi = zeros(Nx,Ny,Np);
    omega = 2*pi*(n-1)/Np*Bandwidth;  % radial frequency
    for ii = 1:Nx
       for jj = 1:Ny
           % zero order phasing
           bestph0 = 0;
           minfom = 1E+6;
%            ph1 = mean(ph1range);
           for ph0=0:.01:2*pi
               spectrum = real(squeeze(complex_img(ii,jj,:)).*exp(1i*(ph0-TE*omega')));
               aux = spectrum;
               aux(aux>0)=0;
               fom=sum(aux.^2);
               if(fom<minfom)
                   minfom=fom;
                   bestph0 = ph0;
%                    mrsi(ii,jj,:) = msbackadj([1:Np]',spectrum);
                    spectrum = real(squeeze(complex_img(ii,jj,:)).*exp(1i*(ph0+pi/10-TE*omega')));%                    mrsi(ii,jj,:) = msbackadj([1:Np]',spectrum);
                    mrsi(ii,jj,:) = msbackadj([1:Np]',spectrum);

               end    
           end
%            if and(ii==9,jj==9)
%            bestph0
%            end
           baseline = polyval(polyfit(baseline_idx,squeeze(mrsi(ii,jj,baseline_idx))',baseline_poly_order),n);
           mrsi(ii,jj,:) = squeeze(mrsi(ii,jj,:))-baseline';
       end
    end
end

%% Show sum spectra
TE = -0.0002;
figure(509)
ph= 5.3400;
sum_spectra = squeeze(sum(sum(mrsi,1),2));
% sum_spectra = real(squeeze(mrsi(8,8,:)).*exp(1i*(ph0-TE*omega')));
baseline = polyval(polyfit(baseline_idx,sum_spectra(baseline_idx)',baseline_poly_order),n);
sum_corrected_spectrum = sum_spectra - baseline';
plot(cs,sum_spectra,cs,sum_corrected_spectrum);
set(gca,'Xdir','reverse');

%% Show MRSI
figure(511)
img_abs = show_csi_matrix_20190409(mrsi,spectrascale,color_peaks,proton,...
    show_proton,xgridshift,ygridshift);

%% Create Maps
% get maps
gas_img = rot90(sum(mrsi(:,:,gas_idx),3),2);
blood_img = rot90(sum(mrsi(:,:,blood_idx),3),2);
tissue_img = rot90(sum(mrsi(:,:,tissue_idx),3),2);
% show maps
figure(512); 
subplot(221);
imagesc(gas_img); colormap hot; axis square; title('gas');
subplot(222);
imagesc(blood_img); colormap hot; axis square; title('blood');
subplot(223);
imagesc(tissue_img); colormap hot; axis square; title('tissue');

% Interpolated maps
gas_imgi = interpolate_image(gas_img,4);
blood_imgi = interpolate_image(blood_img,4);
tissue_imgi = interpolate_image(tissue_img,4);

% Show interpolated maps
figure(513); 
subplot(221);
imagesc(gas_imgi); colormap hot; axis square; title('gas');
subplot(222);
imagesc(blood_imgi); colormap hot; axis square; title('blood');
subplot(223);
imagesc(tissue_imgi); colormap hot; axis square; title('tissue');

% get ratiometric maps
[~,mask] = threshold_image(gas_imgi,0.15);
[~,maskn] = threshold_image(gas_img,0.15);


b2g = blood_img./gas_img.*maskn;
t2g = tissue_img./gas_img.*maskn;
b2t = blood_img./tissue_img.*maskn;


% Chemically-shifted weighted dissolved phase map
for ii = 1:Nx
    for jj =1:Ny
        weighted_s(ii,jj) = sum(squeeze(abs(mrsi(ii,jj,peaks_dis))).*cs_peaks') / sum(squeeze(abs(mrsi(ii,jj,peaks_dis))));    
    end
end

% 
% weighted_gas = fliplr(sum(mrsi(:,:,gas_idx),3).*cs_gas);
% weighted_tissue = fliplr(sum(mrsi(:,:,tissue_idx),3).*cs_tissue);
% weighted_blood = fliplr(sum(mrsi(:,:,blood_idx),3).*cs_blood);

% show ratiometric maps
figure(514); 
subplot(221);
imagesc(b2g);  axis square; title('blood-to-gas'); colormap jet; colorbar;
subplot(222);
imagesc(t2g); colormap jet; axis square; title('tissue-to-gas'); colorbar;
subplot(223);
imagesc(b2t); colormap jet; axis square; title('blood-to-tissue'); colorbar;
colormap jet; caxis([0 1])
subplot(224);
imagesc(fliplr(weighted_s).*maskn); colormap jet; axis square; title('weighted dissolved map'); colorbar;
colormap jet; caxis([185 205])

%%
%% Show KSpace Trajectory
% figure;
% for ii = 1:16
%     for jj = 1:16
%         plot(gx(ii,jj),gy(ii,jj),'o'); hold on;
%         drawnow();
%         pause(0.3); 
%         ylim([0 17])
%         xlim([0 17])
%     end
% end
%% 


