function img = show_csi_matrix_20190409(fnorm7,spectrascale,color_peaks,proton,show_proton,xgridshift,ygridshift,pyr_idx,lac_idx,ala_idx,fov_match,spectral_zoom,xx,yy,xgridscale,ygridscale,xcircshift,ycircshift)
% syntax: img = show_csi_matrix(fnorm,spectrascale,colorpeaks,...
%                               proton,show_proton,...
%                               xgridshift,ygridshift,...
%                               pyr_idx,lac_idx,ala_idx,...
%                               xx,yy,
%                               xgridscale,ygridscale)
%
% This funciton gets the 3D matrix for a 2D csi image and shows it. 
% Given the proton image, it can run overlay the csi image on the proton
% image. 
% Use the output of csirecon(..) as the input of this function
%
% function created by Mehrdad Pourfathi on 2/18/2016.
%
% Update (1) by MP on 2/23/2016
% added a display range in the function itself. 
% changed the color of the spectrum for better visualization.
%
% Update (2) by MP on 3/9/2016
% removed extra variables
% debugged the code.
%
% Update (3) by MP on 6/26/2016
% input arguments "fov_match" and "spectral_zoom" added
% fov_match is set to 1 if carbon and proton fov's are equal.
% spectral_zoom is set to 1 to zoom over the region of the pyr.
% metabolites.
%
% Update (4) by MP on 10/21/2016
% added the "species" input argument. If species is set to 0 then it is
% assumed that the data is from a rat. By default we set the 1H and 13C
% FOVs to 60x60cm2 and 45x45cm2. If species is set to 1, 1H and 13C FOVs
% are assumed to be 40x40cm2 and 25x25cm2.
%
% Update (5) by MP on 2/23/2017
% replaced the "species" input argument with xx and yy that gives that
% actual indicies of the proton image that coincides with the carbon scan
%
% Update (6) by MP on 7/03/2017
% normalizes the input to produce the MRSI map, but does not generate a
% normalized output (I changed it so it does!)
%
% Update (7) by MP on 04/09/2019
% Added "xcircshift" and "ycircshift" as input parameters to be able to
% circularly shift the grid for better visualization as needed.


% manage input parameters
if nargin < 18
    xcircshift = 0;
    ycircshift = 0;
end
if nargin < 16
    xgridscale = 1;
    ygridscale = 1;
end
if nargin < 13
    xx = 1:128;
    yy = 1:128;
end
if nargin < 12
    spectral_zoom = 0;
end
if nargin < 11
    fov_match = 1;
end
if nargin < 8
    pyr_idx = 1;
    lac_idx = 1;
    ala_idx = 1;
end
if nargin < 6
    xgridshift = 0;
    ygridshift = 0;
end
if nargin < 5
    show_proton = 1;
end
if nargin < 4 
    show_proton = 0;
end
if nargin <3
    color_peaks = 0;
end
if nargin < 2
    spectrascale =1;
end

[Nx,Ny,Nt] = size(fnorm7);

hold off;
% show proton image if set to show. 
if show_proton 
    
    imagesc(fliplr((rot90(rot90(proton)))));

    if fov_match 
        [H1Nx,H1Ny] = size(proton);
    else
        [H1Nx,H1Ny] = size(proton(yy(1):yy(end),xx(1):xx(end)));
    end
    
    colormap gray;
    hold on;
else
    H1Nx = Nx;
    H1Ny = Ny; 
%     spectrascale = 1;
end
axis xy          % change axis direction to normal


% grid and spectra color
% colorvec = [0.0431372560560703 0.517647087574005 0.780392169952393];
colorvec = [0 1 1];
% colorvec = [1 1 1];

% setup the variables
img = zeros(Nx,Ny,Nt);

hold on;

% checks if proton and carbon FOVs are matched. 
if fov_match 
    xgrid = xgridscale*linspace(1,H1Nx,Nx+1)+xgridshift;
    ygrid = ygridscale*linspace(1,H1Ny,Ny+1)+ygridshift;
    xline = ygridscale*linspace(1,H1Ny,Ny)+ygridshift;
    yline = xgridscale*linspace(1,H1Nx,Nx)+xgridshift;
else
    xgrid = xx(1)+xgridscale*linspace(1,H1Nx,Nx+1)+xgridshift;
    ygrid = yy(1)+ygridscale*linspace(1,H1Ny,Ny+1)+ygridshift;
    xline = xx(1)+ygridscale*linspace(1,H1Ny,Ny)+ygridshift;
    yline = yy(1)+xgridscale*linspace(1,H1Nx,Nx)+xgridshift;

end


% Create x axis grid
for ii = 1:length(xgrid)
    plot(xgrid(ii)*ones(size(xline)),xline,'Color',colorvec);
end
% Create y axis grid
for ii = 1:length(ygrid)
    plot(yline,ygrid(ii)*ones(size(yline)),'Color',colorvec);
end

% Plot Spectra in the grid;
x = []; y = []; x1 =[]; y1 =[]; x2 = []; y2 = []; x3 = []; y3 = [];
fnorm7 = fnorm7/max(fnorm7(:));             % normalize

fnorm6 = circshift(fnorm7,[xcircshift ycircshift 0]);

fnorm8 = fnorm6;                            % NO! DON'T NORMALIZE
if show_proton
    fnorm8(find(fnorm8(:)>(real(ygrid(2)-ygrid(1))))) = NaN;   % ignore out of range. 
end
spec_cnt = 0;
n = 1:length(fnorm8);

%% create girdline and process spectra
phase_counter = 0;
for ii = 1:Nx
    spec_cnt = spec_cnt + 1;
    xaux = [];
    yaux = [];
    xaux1 =[];
    yaux1 = [];
    xaux2 = [];
    yaux2 = [];
    xaux3 = [];
    yaux3 = [];

    if spectral_zoom
%         disp_range = (lac_idx(1)-10):(pyr_idx(end)+20); % zoom on pyr
        if size(fnorm8,3) == 128
            disp_range = 30:90;
        elseif size(fnorm8,3) == 256
            disp_range = 90:160;
        elseif size(fnorm8,3) == 96
            disp_range = 22:68;
        else
            disp_range = 200:320;
        end
    else
        disp_range = 1:size(fnorm8,3); % show the entire spectrum
%         disp_range = 1:361; % show the entire spectrum

    end
    
    for jj = 1:Ny     
        vec = linspace(xgrid(Nx-jj+1),xgrid(Nx-jj+2),length(disp_range));
%         vec = interp(linspace(xgrid(Nx-jj+1),xgrid(Nx-jj+2),length(disp_range)),4);
        xaux = [xaux NaN vec];
%         xaux = [xaux NaN interp(vec(50):vec(100))];
%         xaux1 = [xaux1 NaN vec(pyr_idx)];
%         xaux2 = [xaux2 NaN vec(lac_idx)];
%         xaux3 = [xaux3 NaN vec(ala_idx)];      
        spec = spectrascale*squeeze(fnorm8(ii,jj,:));
        img(ii,jj,:) = spec*max(fnorm7(:))/spectrascale;
        yaux = [yaux NaN spec(disp_range)'];
%         yaux1 = [yaux1 NaN spec(pyr_idx)']; % show pyruvate peak
%         yaux2 = [yaux2 NaN spec(lac_idx)']; % show lactate peak
%         yaux3 = [yaux3 NaN spec(ala_idx)']; % show alanine peak  
    end 

  
    x = [x NaN xaux];
%     
    yshift = 0;
%     if ii < Ny-yshift
        y = [y NaN yaux+ygrid(ii+yshift)];
%     else
%     	y = [y NaN yaux+ygrid(ii-Ny+yshift+1)];
%     end

%     x1 = [x1 NaN xaux1];
%     y1 = [y1 NaN yaux1+ygrid(ii)];
%     x2 = [x2 NaN xaux2];
%     y2 = [y2 NaN yaux2+ygrid(ii)];
%     x3 = [x3 NaN xaux3];
%     y3 = [y3 NaN yaux3+ygrid(ii)];
end


plot(x,y,'Color',colorvec,'linewidth',2);
axis square;

% if color_peaks
%     plot(x1,y1,'g');
%     plot(x2,y2,'r');
%     plot(x3,y3,'m');
% end
hold off;
