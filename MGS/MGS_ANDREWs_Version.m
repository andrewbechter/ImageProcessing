%test_mgs
clear all
close all

%% load images
%load('NIC_files.mat')
load('NIC2.mat')

%% lbt system parameters
lambda = 1064e-9;

alpha = 0.11;

fnum = 41.2; % system f number

pixel_size = 6.5e-6;

npup = 512; %full size of grid

diameter = 4.4; % just used for f number

focal_length = fnum*diameter; %in mm use this for now

%% create and populate model

model.prior = zeros(size(npup)); % no prior estimate

model.pixel_size = pixel_size; % (pixel size in meters)

model.focal_length = focal_length; % (focal length in mm)

model.diameter = diameter; % (beam diameter in mm)

model.f_number = fnum; % (f number)


%% set processing options
p_opts.rescale = false; %resamples PSF to critical sampling

p_opts.filter = 2; % 1 = fourier filter 2 = weiner

p_opts.psf_rad = 1000; % size of PSF mask [] means no mask

p_opts.diagnostic = false;

p_opts.background = true; 

p_opts.flat = false;

p_opts.cliptol = 1000; % 0 turns off clipping

%% calibrate images
img = calibrate_frame_Andrew(img,p_opts); % fix -ve pixels, divide flats etc...

img = process_psf_Andrew(img,npup,npup,model,p_opts); % filter,mask,centeroid,shift,co-add,resample

%% calculate diversity

if p_opts.rescale
    q = 2; % (assign after resampling option in process_psf)
else
    q = (lambda*model.f_number) * (1/model.pixel_size);
end

% Make a circular pupil (primary)
pup1 = circ(npup,npup/q/2);

% Make a circular pupil (secondary)
pup2 = circ(npup,alpha.*(npup/q/2));

pup = pup1-pup2;

%Calculate amplitude mask after image processing
model.amplitude = pup;

z4 = zernike_mode(pup,4); % Focus mode for diversity

%dz = extractfield(img,'score_z'); % get stage positions in microns

dz = [0,12,6,-6,-12];

dz = dz.*1e3*1e-6; %adjust LBT stage positions to be in microns

div_PV = calc_div(dz,lambda,fnum); % calculate diversity (P-V?) in waves (dz and lambda in same units)

divamp = 0.5*div_PV; %amplitude is half the P-V for focus...

for k = 1:length(divamp)
    
    img(k).diversity = divamp(k)*lambda.*(z4./(max(max(z4)))); % diversity amplitudes in meters
   
end

temp = img;

samp_psf = img(1).psf;

img(1)=[]; % remove in focus psf

%% set mgs options
opts.adapt_diversity = false;
opts.adapt_tilt = true;
opts.adapt_iter = 2;
opts.adapt_method = 'mgs';
adapt_options=[];
adapt_zernikes=[4,5,6];
adapt_tol = [1e-9];

opts.num_prior_loops=1;
opts.num_outer_loops=7;
opts.num_inner_loops = 200;
opts.convergence_tol = 1e-20;

[wfs,adapt] = mgs(img,model,opts);

%% run non-adaptive MGS
% wfs = mgs_matlab(img,model,opts);

wfs.mgs.opd = zernike_remove(wfs.mgs.opd,wfs.mgs.opd~=0,1); %remove

fprintf('Estimated OPD: %g nm\n',std(nonzeros(wfs.mgs.opd))*1e9);
fprintf('Estimated Strehl: %g\n',(exp(-(2*pi*(std(nonzeros(wfs.mgs.opd)))/lambda)^2)));

%% figures

rec_psf = opd2psf(pup,wfs.mgs.opd,lambda);

fprintf('O-C: %g \n',std(nonzeros(rec_psf-samp_psf)));

pos_num =1;

psf = img(pos_num).psf; % assign the measured and calibrated in focus image to the psf

mgs_psf = wfs.mgs.psf_est(:,:,pos_num);

figure
subplot(2,3,1)
imagesc(pad(psf,150,150))
axis image
title('PSF')
colorbar
cl = caxis;
axis image

subplot(2,3,2)
imagesc(pad(mgs_psf,150,150))
title('MGS PSF')
colorbar
caxis(cl)
axis image

subplot(2,3,3)
imagesc(pad(wfs.mgs.opd,150,150));
title('OPD estimate');
axis image
colorbar;
caxis([-6e-7,6e-7])

subplot(2,3,4)
imagesc(pad(log(abs(samp_psf)),150,150))
title('Real "in focus" PSF sample frame')
axis image
colorbar;

subplot(2,3,5)
imagesc(pad(log(abs(rec_psf)),150,150))
axis image
title('OPD to PSF')
colorbar;

subplot(2,3,6)
imagesc(pad(log(abs(rec_psf-samp_psf)),150,150))
title('PSF O-C')
colorbar
axis image

function img = calibrate_frame_Andrew(img,p_opts)
% this function will eventually do all the calibration if necessary (right
% now it only fixes -ve pixels and is kind of pointless)

if nargin<2
    p_opts.diagnostic = false;
    p_opts.background = false;
end

for n = 1:length(img)
    
    frames  = img(n).frames; %assign frames at one location
    
    for f = 1:size(frames,3)
        
        p = frames(:,:,f);
        
        %% -------------------------------------
        % Remove background (aperture method)
        % --------------------------------------
        if p_opts.background == true
            
            if p_opts.diagnostic == true && f==1
                figure(1);
                imagesc(log(abs(p)))
                title('raw')
                colorbar
                pause(1)
            end
            
            [b,~] = background(p);
            
            p = p-b;
            
            if p_opts.diagnostic==true && f==1
                figure(1);
                imagesc(log(abs(p)))
                title('background')
                colorbar
                pause(1)
            end
            
        end
        
        if p_opts.flat == true
           if isempty(p_opts.flatname)
               p_opts.flatname = '/Users/bechter/Documents/MATLAB/Andrew_MGS_060818/masterFlat_20180819b.mat';
           end
               
            flat = load(p_opts.flatname);
            flat = flat.master_Flat;
            
            p = p./flat;
            
            if p_opts.diagnostic == true && f==1
                figure(1);
                imagesc(log(abs(p)))
                title('flat')
                colorbar
                pause(1)
            end
        end
            
        
%         %% ---------------------------
%         % Clip Negative pixel values 
%         % ----------------------------
%         
%         [ind]= find(p<0);
%         
%         p(ind)=0;
%         
         frames(:,:,f) = p;
         
    end
    
    img(n).frames = frames;

end
end
function [b,bsig] = background(frame)

imageSizeX = size(frame,2); % size of frame cols

imageSizeY = size(frame,1); % size of frame rows

[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);

% Next create the apreture circle in the image.

centerX = imageSizeX/2;

centerY = imageSizeY/2;

innerRadius = imageSizeX/4;

outerRadius = imageSizeX/3;

array2D = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2;

ringPixels = array2D >= innerRadius.^2 & array2D <= outerRadius.^2;

% ringPixels is a 2D "logical" array.

b = mean(frame(ringPixels));

bsig = std(frame(ringPixels));

end
function out = pad(in,npix_rows,npix_cols)
% PAD zero-pads input matrix to size npix by npix or cuts down to size\
%
% out = pad_truncate(in,npix_rows,npix_cols)
%
% Parameters
% ----------
% in : m x n array
%   Matrix to be padded
%
% npix_rows : int
%   Desired number of rows in output
% 
% npix_cols : int, optional
%   Desired number of columns in output. Default is npix_rows (result will
%   be square: npix_rows x npix_rows)
%
% Returns
% -------
% out : npix_rows x npix_cols array
%   Padded input 
%

narginchk(2,3)

if nargin < 3
    npix_cols = npix_rows;
end

out = zeros(npix_rows,npix_cols);

[nrows, ncols]  = size(in);

ixc = floor(ncols/2 + 1);
iyc = floor(nrows/2 + 1);

oxc = floor(npix_cols/2 + 1);
oyc = floor(npix_rows/2 + 1);

dx = npix_cols-ncols;
dy = npix_rows-nrows;

if dx<=0
    ix1 = ixc - floor(npix_cols/2);
    ix2 = ix1 + npix_cols - 1;
    ox1 = 1;
    ox2 = npix_cols;
else
    ix1 = 1;
    ix2 = ncols;
    ox1 = oxc - floor(ncols/2);
    ox2 = ox1 + ncols - 1;
end

if dy<=0
    iy1 = iyc - floor(npix_rows/2);
    iy2 = iy1 + npix_rows - 1;
    oy1 = 1;
    oy2 = npix_rows;
else
    iy1 = 1;
    iy2 = nrows;
    oy1 = oyc - floor(nrows/2);
    oy2 = oy1 + nrows - 1;
end

out(oy1:oy2, ox1:ox2)  = in(iy1:iy2, ix1:ix2);
end
function [img_centered, dx, dy] = centroid_center(img, threshold, niter, verbose)
% center iteratively centers a PSF image using centroids
%
% [img_centered, dx, dy] = CENTROID_CENTER(img, threshold, niter, verbose*)
%
% Parameters
% ----------
% img : m x n array
%   PSF image
%
% threshold : float
%   Clipping threshold relative to peak
%
% niter : int
%   Maximum number of iterations to use
%
% verbose : bool, optional
%   If true, print details to the console. Default is false.
%
% Returns
% -------
% img_centered : m x n array
%   Centered and unit normalized PSF
%
% dx : float
%   Measured x-offset (column offset)
%
% dy : float
%   Measured y-offset (row offset)
%

narginchk(3,4)
if nargin < 4
    verbose = false;
end

[nrows, ncols] = size(img);

xc = floor(ncols/2) + 1;
yc = floor(nrows/2) + 1;

img_centered = img;

dxh=[];
dyh=[];

for kk = 1:niter
    [xo, yo] = centroid(img_centered.*(img_centered>max(img_centered(:)*threshold)));
    
    dx = xo-xc;
    dy = yo-yc;
    
    dxh = [dxh dx];
    dyh = [dyh dy];
    
    img_centered = image_shift(img_centered, -dx, -dy);
    
    dr = sqrt(dx.^2 + dy.^2);
    
    vprintf(verbose, ' Iter=%d : dx=%5.2f dy=%5.2f (delta_iter=%6.3f)\n',kk,sum(dxh),sum(dyh),dr)
    
    if dr<0.001
        break;
    end
end

dx = sum(dxh);
dy = sum(dyh);
vprintf(verbose, ' PSF centered to within %6.4f pixels.\n',dr)
vprintf(verbose, '------------------------------------------------------------\n ')

end
function [img] = process_psf_Andrew(img,nrow,ncol,model,p_opts)
% break out model and img constants and assign function specific parameters
threshold = 0.5;

iterations = 5;

snr_tol = p_opts.cliptol;

verbose = false;

lambda = img(1).bandpass(1);

pixel_size = model.pixel_size;

f_number = model.f_number;

q = (lambda*f_number) * (1/pixel_size);

%% Trim/pad, then centroid and shift each frame to center on PSF (coadds frames at each z)

for n = 1:length(img)
    
    psum = 0;
    
    frames  = img(n).frames; %assign frames at one location
    
    for f = 1:size(frames,3)
        
        p  = frames(:,:,f); %assign indiviudal frames
        
        %% -------------------------------------
        % Pad to correct  rectangular frame
        % --------------------------------------
        
        dim1 = min(size(p));
        
        p = pad(p,dim1,dim1); % pad single frame
        
        
        %% -----------------------------------------
        % Filter (fourier method = 1) (weiner2d = 2)
        % -----------------------------------------
        
        if p_opts.filter ==1
            %q = fraction of frame to mask
            q = 10; %testing this out now
            
            p = f_filter(p,q);
            
        elseif p_opts.filter==2
            
            p = wiener2(p,[2,2]);
            
        end
        
        if p_opts.diagnostic && f==1
            figure(1);
            imagesc(log(abs(p)))
            title('filter')
            colorbar
            pause(1)
        end
        
        %% -------------------------------------
        % Apply large PSF mask to remove ghosts(circle method)
        % --------------------------------------
        npix = size(p,1);
        
        if ~isempty(p_opts.psf_rad)
            xs = [-(npix/2):(npix/2-1)];
            [x,y] = meshgrid(xs,xs);
            psf_mask = sqrt(x.^2 + y.^2) < p_opts.psf_rad;
        else
            psf_mask = ones(npix);
        end
        
        p = p.*psf_mask;
        
        if p_opts.diagnostic && f==1
            figure(1);
            imagesc(log(abs(p)))
            title('Ghost mask')
            colorbar
            pause(1)
        end
        
        %% -------------------------------------
        % Center frames
        % --------------------------------------
        
        p = centroid_center(p, threshold, iterations,verbose);
        
        if p_opts.diagnostic && f==1
            figure(1);
            imagesc(log(abs(p)))
            title('center')
            colorbar
            pause(1)
        end
        
        %% -----------------------------------------
        % Trim/Pad
        % -----------------------------------------
        
        p = pad(p,nrow,ncol); % pad single frame
        
        if p_opts.diagnostic && f==1
            figure(1);
            imagesc(log(abs(p)))
            title('pad')
            colorbar
            pause(1)
            
        end
        
        npix = size(p,1);
        
        %% -------------------------------------
        % Center frames
        % --------------------------------------
        
        p = centroid_center(p, threshold, iterations,verbose);
        
        if p_opts.diagnostic && f==1
            figure(1);
            imagesc(log(abs(p)))
            title('center')
            colorbar
            pause(1)
        end
        
        % -----------------------------------------
        % Co-add frames
        % -----------------------------------------
        
        psum = psum + p;
        
    end
    
    if p_opts.diagnostic
        figure(1);
        imagesc(log(abs(psum)))
        title('sum')
        colorbar
        pause(1)
    end
    %% ------------------------------------
    % Resample PSF to be critically sampled
    % -------------------------------------
    if p_opts.rescale % rescales to q = 2
        
        psum = rescale_psf(psum,lambda,f_number, pixel_size, npix);
   
    psum = pad(psum,npix);
    
    if p_opts.diagnostic 
        figure(1);
        imagesc(log(abs(psum)))
        colorbar
        pause(1)
    end
    
     end
    
    % -----------------------------------------
    % Fix negative pixels before normalizing
    % -----------------------------------------
    psum(psum<0) = 0;
    
    if p_opts.diagnostic
        figure(1);
        imagesc(log(abs(psum)))
        title('Clip negatives')
        colorbar
        pause(1)
    end
    
  
    
    
    %% -------------------------------------
    % Apply small PSF mask to remove noise
    % --------------------------------------
    npix = size(p,1);
    
    psf_rad = [];
    
    if ~isempty(psf_rad)
        xs = [-(npix/2):(npix/2-1)];
        [x,y] = meshgrid(xs,xs);
        psf_mask = sqrt(x.^2 + y.^2) < psf_rad;
    else
        psf_mask = ones(npix);
    end
    
    psum = psum.*psf_mask;
    
    if p_opts.diagnostic
        figure(1);
        imagesc(log(abs(psum)))
        colorbar
        title('Noise-mask')
        pause(1)
    end
    
    %% -----------------------------------------
    % Clip PSF intensities below SNR tolerance
    % ------------------------------------------
    
    if snr_tol ~= 0
        snr_threshold = max(max(sqrt(psum.*(psum>0)))) / snr_tol;
        snr_mask = sqrt(psum.*(psum>0)) > snr_threshold;
        psum = psum .* snr_mask;
        
        
        if p_opts.diagnostic
            figure(1);
            imagesc(log(abs(psum)))
            colorbar
            title('Clip low SNR pixels')
            pause(1)
        end
        
    end
    
    disp(sum(sum(psum)))
    % -----------------------------------------
    % Normalize intensity
    % -----------------------------------------
    psum = psum./sum(psum(:));
    
    if p_opts.diagnostic
        figure(1);
        imagesc(log(abs(psum)))
        colorbar
        title('Normalize')
        pause(1)
    end
    
    img(n).psf = psum;
end

end
function [xc, yc] = centroid(img)
% CENTROID computes image centroid location
%
% [xc, yc] = CENTROID(img)
%
% Parameters
% ----------
% img : m x n array
%   Image
%
% Returns
% -------
% x : float
%   x centroid (column)
%
% y : float
%   y centroid (row)
%

[rows, cols] = size(img);

xc = floor(cols/2) + 1;
yc = floor(rows/2) + 1;

m  = sum(sum(img));

if m ~= 0
    mx = sum(img)  * (1:cols)';
    my = sum(img') * (1:rows)';

    xc = mx / m;
    yc = my / m;
end
end
function [psf_new] = rescale_psf(psf_image, lambda, f_num, pix_size, winsize)

lambda_nyq = 2 * pix_size / f_num;
npix       = max(size(psf_image));
npix2      = round(npix * lambda_nyq/lambda);

g1 = [1:npix]./npix;
g2 = [1:npix2]./npix2;
[x1 y1] = meshgrid(g1,g1);
[x2 y2] = meshgrid(g2,g2);

psf_new = pad(interp2(x1,y1,pad(psf_image,npix),x2,y2),winsize);

psf_new(find(isnan(psf_new))) = 0;

end
function [wfs,adapt] = mgs(img,model,options)

%% ---------------------------------------
% Parse inputs and check structure fields
% ----------------------------------------

if nargin==2
    options = mgs_options;
end

%% -----------------------------------------------------------------
% Preprocess images 
% Note: for now, pupil and diversity OPDs generated outside wfsc-mgs
% -------------------------------------------------------------------

if ~isfield(img,'psf')
    error('No PSF images in IMG structure! See IMG data structure definition.');
end


%% ---------------------------------------
% Diversity adaptation
% ----------------------------------------

if options.adapt_diversity
    
    % Run adaptation
    adapt = mgs_adapt_diversity(img,model,options);

    % Add adaptation phase to diversity
    for k=1:length(img)
        img(k).diversity = img(k).diversity +  adapt(k).tt + adapt(k).opd;
    end
    
else
    adapt = [];
end

%% ---------------------------------------
% Run MGSs
% ----------------------------------------

m = mgs_matlab(img,model,options);


%% ---------------------------------------
% Populate output structure
% ----------------------------------------

wfs.mgs = m;
wfs.img = img;
wfs.options = options;
end
function adapt = mgs_adapt_diversity(img,model,options)
% div = mgs_adapt_diversity(img,model,options)
% Runs MGS on each PSF in IMG to estimate the non-common phase.
%
%
% Inputs:
% -------
% img : struct
%   Standard WFSC image structure. See wfsc-dev data structures definition.
%
% model : mgs_model object
%   MGS model object. This can be instantiated with mgs_model.
%
% options : mgs_options object
%   MGS options object. This can be instantiated with default values with
%   mgs_options.
%
%
% Outputs:
% --------
% adapt : struct
%   Structure containing the adaptation phase for each PSF in img
%   adapt(k).tt     Tip/tilt OPD for img(k)
%   adapt(k).opd    Adaptation OPD for img(k)

%% -------------------------------------
% Extract values from input structures
% --------------------------------------

adapt_tilt = options.adapt_tilt;
adapt_iter = options.adapt_iter;
adapt_tol = options.adapt_tol;
adapt_method = options.adapt_method;
zmodes = options.adapt_zernikes;

if ~isfield(options,'est_opts') || isempty(options.est_opts)
    warning('Using default MGS options for diversity adaptation');
    mgs_opts = options;
else
    mgs_opts = options.adapt_options;
end

mask = model.amplitude~=0;
f_number = model.f_number;
px = model.pixel_size;

[gr,ga] = zernike_coordinates(mask,'CIRCLE');

%% --------------------------------------
% Global adaptation
% ---------------------------------------

% Setup tilt modes (unit RMS)
z2 = zernike_mode(mask,2,gr,ga);
z3 = zernike_mode(mask,3,gr,ga);

for n=1:length(img)

    imgn = img(n);              % this image

    % Estimate global tilt from PSF centroid
    % 8 Nyq Pix = 1 lambda rms
    if adapt_tilt
        p = imgn.psf;
        lambda = imgn.bandpass(1,1);
        q = lambda*f_number/px;
        [~,xc,yc] = centroid_center(p,0.005,50);
        zk2 = (2/q)*(xc*lambda/8);
        zk3 = (2/q)*(yc*lambda/8);
        tt = zk2*z2 + zk3*z3;
    else
        tt = 0;
    end
    imgn.diversity = imgn.diversity + tt;   % Add tt term to diversity

    % Run adaptation iterations
    for m=1:adapt_iter
        switch lower(adapt_method)
            case 'mgs'
                wfs = mgs_matlab(imgn,model,mgs_opts);
            case 'ipo'
                error('IPO adaptation not yet supported!');
            otherwise
                error('Unrecognized adaptation method!');
        end
        [fit,zk] = zernike_fit(wfs.opd, mask, zmodes, gr, ga);        
        imgn.diversity = imgn.diversity + fit;
        tol_flag = ~isempty(find(abs(zk') > adapt_tol)); 
        if tol_flag == 0 break; end
    end

    adapt(n).tt = tt;                                       % TT adapt OPD
    adapt(n).opd = imgn.diversity - tt - img(n).diversity;  % non-TT adapt OPD
end
end
function [rho,theta] = zernike_coordinates(mask,method,center)
% ZERNIKE_COORDINATES generates a Zernike coordinate system
%
% [rho,theta] = ZERNIKE_COORDINATES(mask,method,center)
%
% Parameters
% ----------
% mask : m x n array
%   Binary mask defining the extent to compute the Zernike polynomials over.
%
% method : str, optional
%   Coordinate normalization method. Default is 'CIRCLE'. Valid values are:
%
%       'CIRCLE'  - outscribing circle radius =  1
%                   center coordinates = 0
%       'AXIS'    - largest X departure = 1
%                   largest Y departure = 1
%                   center coordinates = 0
%       'GENERAL' - exterior perimeter = 1
%                   interior perimeter = 0
%
% center : 1 x 2 array, optional
%   [x,y] center in pixels. Default is the mask centroid
%
% Returns
% -------
% rho : m x n array
%   Radial coordinates of the mask array.
%
% theta : m x n array
%   Angular coordinates of the mask array in radians.
%


narginchk(1,3)

mask = (mask~=0);
[rows, cols] = size(mask);

if nargin < 2
    method = 'CIRCLE';
end

% Determine the Central Ordinate of the Pupil
if nargin < 3
    [xc, yc] = centroid(mask);
else
    xc = center(1);
    yc = center(2);
end

% Define Initial Pixel Coordiate System
xs = (1:cols)-xc;
ys = (1:rows)-yc;
[x, y] = meshgrid (xs,ys);

r = abs(x + 1j*y);
a = angle(x + 1j*y);

rmin = 0;
rmax = max(r(:).*mask(:));

switch upper(method)
    case 'CIRCLE'
        rho   = r ./ rmax;
        theta = a;
        
    case 'AXIS'
        dx = max(abs(x(:).*mask(:)));
        dy = max(abs(y(:).*mask(:)));
        rho   = abs(x./dx + 1j*y./dy);
        theta = angle(x./dx + 1j*y./dy);
        
    case 'GENERAL'
        rmax = 0*r;
        rmin = 0*r;
        
        % Determine The rough boundaries of the pupil
        
        probe = (-pi:pi/rows*4:pi-pi/rows*4);
        pmin = [];
        pmax = [];
        mm = 0;
        for n=1:length(probe)
            amask = abs((a+probe(n)).*r)<1;
            rlist = find (amask.*mask);
            pmin(n) = min(r(rlist));
            pmax(n) = max(r(rlist));
            mm = mm + (amask.*mask);
        end
        
        probe = [(probe-2*pi) probe (probe+2*pi)];
        pmax  = [pmax pmax pmax];
        pmin  = [pmin pmin pmin];
        
        list = find(mask);
        for n=1:length(list)
            ai = a(list(n));
            ri = r(list(n));
            
            rmax(list(n)) = interp1(probe,pmax,ai);
            rmin(list(n)) = interp1(probe,pmin,ai);
        end
        
        rho   = (r-rmin) ./ (rmax-rmin);
        rho(find(isnan(rho)))=0;
        %rho = rho.* (rho<=1) + (rho>1);
        theta =  a;
        
    otherwise
        fprintf('zernike_coordinates: Method = %s not recognized.\n', method);
        
end
end
function [zern_opds,m,n] = zernike_mode(mask,modes,rho,theta)
% ZERNIKE_MODE computes the circular Noll Zernike polynomial for a given mask.
%
% [zern,m,n] = zernike_mode(mask,modes,rho,theta)
%
% Parameters
% ----------
% mask : m x n array
%   Binary mask defining the extent to compute the Zernike polynomials over.
%
% modes : int or 1 x N array
%   Noll ordered mode(s) to create.
%
% rho : m x n array, optional
%   Radial coordinates of the mask array. rho should be 0 at the origin and 1
%   at the edge of the circle. Default is computed for an outscribing circle
%   over the extent of the mask.
%
% theta : m x n array, optional
%   Angular coordinates of the mask array in radians. Default is computed for an
%   outscribing circle over the extent of the mask.
%
% Returns
% -------
% zern : m x n x N array
%   Zernike mode(s)
%
% m, n : int
%   Azimuthal and radial indices
%


if nargin < 4
    [rho, theta] = zernike_coordinates(mask, 'CIRCLE');
end

% ---------------------------------------------------
% Compute m and n for the particular zernike mode
% ---------------------------------------------------

j = 1;
zm = [];
zn = [];
for n = 0:round(max(modes)/2)
    for m = (0:2:n)+(mod(n,2)~=0)
        if m > 0
            zm(j) = m;
            zn(j) = n;
            j = j + 1;
            zm(j) = m;
            zn(j) = n;
            j = j + 1;
        else
            zm(j) = m;
            zn(j) = n;
            j = j + 1;
        end
    end
end

% ---------------------------------------------------
% Compute zernike mode over pupil support
% ---------------------------------------------------

iloop = 0;

for imode = modes
    
    m = zm(imode);
    n = zn(imode);
    
    R = 0*mask;
    for s=0:((n-m)/2)
        Rk = (-1)^s*factorial(n-s)/factorial(s)/factorial((n+m)/2-s)/factorial((n-m)/2-s);
        R = R + Rk * rho.^(n-2*s);
    end
    
    zern = sqrt(n+1) * R;
    
    if m > 0
        if mod(imode,2)==0
            zern = zern .* cos(m*theta) .* sqrt(2);
        else
            zern = zern .* sin(m*theta) .* sqrt(2);
        end
    end
    
    iloop = iloop + 1;
    zern_opds(:,:,iloop) = zern .* (mask~=0);
    
end
end
function [opd_fit,zern_coeff] = zernike_fit(opd,mask,modes,rho,theta)
% ZERNIKE_FIT computes Zernike coefficients of a given OPD
%
% [opd_fit,zern_coeff] = ZERNIKE_FIT(opd,mask,order,rho,theta)
%
% Parameters
% ----------
% opd : m x n array
%   OPD to be fit
%
% mask : m x n array 
%   Binary mask defining the extent to fit the Zernike polynomials over
% 
% modes : int or 1 x N array
%   Noll ordered mode(s) to be fit
%
% rho : m x n array, optional
%   Radial coordinates of the mask array. rho should be 0 at the origin and 1 
%   at the edge of the circle. Default is computed for an outscribing circle 
%   over the extent of the mask.
%
% theta : m x n array, optional
%   Angular coordinates of the mask array in radians. Default is computed for an 
%   outscribing circle over the extent of the mask.
%
% Returns
% -------
% opd_fit : m x n array
%   OPD composed of the fitted terms
%
% zern_coeff : float or 1 x N array
%   Noll ordered Zernike coefficients fit to opd
%

narginchk(3,5)
if nargin == 3
    zern = zernike_mode(mask,modes);
elseif nargin == 5
    zern = zernike_mode(mask,modes,rho,theta);
else
    error('You must specify both rho and theta')
end

zern_coeff(modes) = mode_fit(opd,mask,zern);
opd_fit = zernike_compose(mask,zern_coeff);
end
function opd = zernike_compose(mask,zern_coef,rho,theta)
% ZERNIKE_COMPOSE creates an OPD based on Zernike coefficients.
%
% opd = ZERNIKE_COMPOSE(mask,zern_coef,rho,theta)
%
% Parameters
% ----------
% mask : m x n array
%   Binary mask defining the extent to compute the Zernike polynomials over.
%
% zern_coef : 1 x n array
%   Coefficients corresponding to Zernike indices (Noll ordering) used to
%   create the OPD.
%
% rho : m x n array, optional
%   Radial coordinates of the mask array. rho should be 0 at the origin and 1 
%   at the edge of the circle. Default is computed for an outscribing circle 
%   over the extent of the mask.
%
% theta : m x n array, optional
%   Angular coordinates of the mask array in radians. Default is computed for an 
%   outscribing circle over the extent of the mask.
%
% Returns
% -------
% opd : m x n array
%
% Notes
% -----
% Coefficients less than 1e-12 are ignored in the composition. 
%

opd = zeros(size(mask));
zi = find(abs(zern_coef) > 1e-12);

if isempty(zi)
    return
end

for n=1:length(zi)
    if nargin < 3
        opd = opd + zern_coef(zi(n)) * zernike_mode(mask, zi(n));
    else
        opd = opd + zern_coef(zi(n)) * zernike_mode(mask, zi(n), rho, theta);
    end
end
end
function coeff = mode_fit(opd,mask,modes)
% MODE_FIT fits modes to an OPD by solving the system of linear equations
% modes * coeff = opd for coeff.
%
% coeff = MODE_FIT(opd,mask,modes)
%
% Parameters
% ----------
% opd : m x n array
%   OPD to be fit
%
% mask : m x n array 
%   Binary mask defining the extent to fit supplied modes over
% 
% modes : m x n x N array
%   Cube of mode shapes to fit
%
% Returns
% -------
% coeff : N x 1 array
%   Coefficients fit to OPD
%


% vectorize everything
mask = mask(:);
opd = opd(:);
[mr, mc, mn] = size(modes);
modes = reshape(modes,mr*mc,mn);

% apply mask
opd = opd(mask ~= 0);
modes = modes(mask ~= 0,:);

% fit modes to opd
coeff = modes\opd;
end
function wfs = mgs_matlab(img,model,options)
persistent cvec

% MGS_MATLAB Runs MGS in Matlab for the images contained in the img input structure.
%
% wfs = MGS_MATLAB(img,pup_model,mgs_options) runs the Modified Gerchberg-Saxton phase % retrieval algorithm on the PSF images contained in img using the specified
% model and options structures. This version of MGS is entirely Matlab-based.
%
% Inputs:
% -------
% img : struct
%   Standard WFSC image structure.
%
% pup_model : struct
%   Structure containing the following fields:
%   .prior      Prior OPD estimate
%   .amplitude  Pupil amplitude map
%
% mgs_options : struct
%   Structure containing the following fields:
%   .num_inner_loops    Number of MGS loops for each diversity PSF
%   .num_outer_loops    Number of MGS loops over set of diversity PSFs
%   .num_prior_loops    Number of iterations with prior updating
%
%
% Outputs:
% --------
% wfs : struct
%   Structure containing the phase retrieval result.
%   .opd_mgs    Differential (from prior) pupil phase estimate.
%   .opd        Total pupil phase estimate (opd_mgs + prior)
%   .psf_est    Matched PSF images at each diversity location
%   .cost       Cost function history
%   .diff       Difference OPD history
%
%
% References:
% -----------
% [1] Bikkannavar, S., et al., "Phase Retrieval Methods for Wavefront Sensing,"
% Proc. of SPIE, 2010.
%
% [2] Green, J., et al., "Extreme Wave Front Sensing Accuracy for the Eclipse
% Coronographic Space Telescope," Proc. of SPIE, 2003.

% Jonathan Tesch, Jet Propulsion Laboratory
% Copyright 2018 California Institute of Technology

tic;


%% -------------------------------------
% Extract values from input structures
% --------------------------------------

% img
if size(img(1).bandpass,1) > 1
    warning('Multiwavelength MGS not supported, using maxium wavelength');
    [~,idx] = max(img(1).bandpass(:,2));
    lambda = img(1).bandpass(idx,1);
else
    lambda = img(1).bandpass(1,1);
end
num_psf = length(img);
num_div = length(img);
npix_psf = max(size(img(1).psf));
npix_div = max(size(img(1).diversity));

% Model
prior0 = model.prior;
pupil = model.amplitude;
mask = model.amplitude~=0;
npix_pup = max(size(pupil));
npix_opd = max(npix_pup,npix_div);
if npix_pup > npix_opd
    pupil = pad(pupil,npix_opd);
    mask = pupil~=0;
    npix_pup = max(size(pupil));
end

% MGS options
nInnerLoops = options.num_inner_loops;
nOuterLoops = options.num_outer_loops;
nPriorLoops = options.num_prior_loops;
convTol = options.convergence_tol;
if ~isfield(options,'phase_limit')
    qlimit = pi*ones(npix_opd);
else
    qlimit = options.phase_limit*ones(npix_opd);
end


%% -------------------------------------
% Initialize images
% --------------------------------------

% Pad diveristy phases
div = zeros(npix_opd,npix_opd,num_div);
for k=1:num_div
    div(:,:,k) = pad(img(k).diversity,npix_opd);
end

% Normalize images
psf = zeros(npix_psf,npix_psf,num_psf);
for k=1:num_psf
    psf(:,:,k) = img(k).psf/sum(img(k).psf(:));
end

% Initial phase estimate
rphase_m = zeros(npix_opd,npix_opd);		% joint estimate over all diversity phases



%% -------------------------------------
% Run MGS
% --------------------------------------

fprintf('----------------------------------------------\n');
fprintf('Running MATLAB MGS\n')
fprintf('Number of diversity OPDs: %d\n',num_div);
fprintf('Prior loops: %d\n',nPriorLoops);
fprintf('Outer loops: %d\n',nOuterLoops);
fprintf('Inner loops: %d\n',nInnerLoops);
fprintf('----------------------------------------------\n');

cvec = [];		% accumulate cost function
evec = [];		% accumulate difference from prior
prior = prior0;
for kk = 1:nPriorLoops
    
    for jj = 1:nOuterLoops
        fprintf('Outer loop %d of %d\n',jj,nOuterLoops);
        
        % Loop over each diversity OPD
        rphase = zeros(npix_opd,npix_opd,num_div);
        psf_est = zeros(npix_psf,npix_psf,num_div);
        for ii = 1:num_div
            qest = rphase_m;
            qdiv = div(:,:,ii);
            qpri = prior;
            
            for iter = 1:nInnerLoops
                pfield0 = pupil .* exp(2*pi/lambda*1j*(qest+qdiv+qpri));	% pupil field
                ifield0 = ifftshift(fft2(pfield0)/npix_psf);				% image field
               
%                 keyboard
%                 figure(100)
%                 subplot(1,3,1)
%                 imagesc(ifield0.*conj(ifield0))
%                 subplot(1,3,2)
%                 imagesc(img(ii).psf)
%                 subplot(1,3,3)
%                 imagesc(img(ii).psf-(ifield0.*conj(ifield0))./sum(sum(ifield0.*conj(ifield0))))

                iphase = angle(ifield0);
                
                ifield1 = sqrt(psf(:,:,ii)) .* exp(1j*iphase);
                pfield1 = ifft2(ifftshift(ifield1)*npix_psf);
                
                % Solve for qest by subtracting off diveristy and prior and thresholding
                qest = angle(pfield1 .* exp(-1j*2*pi/lambda*(qdiv+qpri)))/2/pi*lambda;
                qest = min(qlimit,qest);
                qest = max(-qlimit,qest);
                
                
            end
            
%             figure(11)
%             subplot(1,2,1)
%             imagesc(ifield0.*conj(ifield0))
%             subplot(1,2,2)
%             imagesc(qest)
%             drawnow
%             pause()
            
            rphase(:,:,ii) = qest;
            psf_est(:,:,ii) = conj(ifield0).*ifield0;
            psf_est(:,:,ii) = psf_est(:,:,ii)/sum(sum(psf_est(:,:,ii)));
        end
        
        % Calculate joint estimate
        rphase_m = mean(rphase,3);
        estimate = (rphase_m + prior).*mask;
        
        % Calculate cost functiion: mean square of PSF difference
        diff = prior - estimate;
        diff = diff - mask.*mean(nonzeros(diff));
        
        for k=1:num_div
            pc(k) = sum(sum((psf(:,:,k) - psf_est(:,:,k)).^2));
        end
        psf_cost = mean(pc);
        
        cvec  = [cvec psf_cost];
        
        figure(10)
        plot(cvec)
        drawnow
        
        evec = [evec std(nonzeros(diff))];
        
        fprintf('\tCost: %6.5g\n\n',psf_cost);
        
        %fprintf('\tDiff: %6.5g\n\n',std(nonzeros(diff)));
        
        if psf_cost < convTol
            fprintf('Convergence criteria reached, exiting\n');
            break;
        else
            prior = estimate;
        end
    end
end

% Setup outputs
estimate = pad(estimate,npix_opd);
wfs.opd = estimate + prior0;		% total estimated OPD
wfs.opd_mgs = estimate;				% OPD estimated by MGS
wfs.psf_est = psf_est;				% matched PSFs
wfs.cost = cvec;					% cost evolution
wfs.diff = evec;					% OPD error evolution

toc;
end
function [result,zern_coeff] = zernike_remove(opd,mask,modes,rho,theta)
% ZERNIKE_REMOVE computes and removes Zernike coefficients of a given OPD
%
% [result,zern_coeff] = ZERNIKE_REMOVE(opd,mask,order,rho,theta)
%
% Parameters
% ----------
% opd : m x n array
%   OPD to be fit
%
% mask : m x n array 
%   Binary mask defining the extent to fit the Zernike polynomials over.
% 
% modes : int or 1 x N array
%   Noll ordered mode(s) to be removed
%
% rho : m x n array, optional
%   Radial coordinates of the mask array. rho should be 0 at the origin and 1 
%   at the edge of the circle. Default is computed for an outscribing circle 
%   over the extent of the mask.
%
% theta : m x n array, optional
%   Angular coordinates of the mask array in radians. Default is computed for an 
%   outscribing circle over the extent of the mask.
%
% Returns
% -------
% result : m x n array
%   Residual OPD after removal of the fit terms
%
% zern_coeff : float or 1 x N array
%   Noll ordered Zernike coefficients removed from OPD
%


if nargin == 5
    [zfit, zern_coeff] = zernike_fit(opd,mask,modes,rho,theta);
else
    [zfit, zern_coeff] = zernike_fit(opd,mask,modes);
end

result = opd - zfit;
end
function [D] = calc_div(dz,lambda,fnum)
%[D] = calc_div(dz,lambda,fnum) 
% sets mgs options, runs mgs_matlab and hadles plotting
%
% Inputs:
% -------
% dz : 1xn double
%   vector of stage positions with 0 being nominal focus (in meters)
%
% lambda : double
%   wavelength in meters (or same units as dz)
%
% fnum : double
%   score f number
%
% Outputs:
% --------
% D : peak to valley diversity in waves

D = dz./(8 * lambda * (fnum)^2); % only for focus mode 

end
function result = image_shift(image,dx,dy,complex)
% IMAGE_SHIFT shifts an image via FFT 
%
% result = IMAGE_SHIFT(image,dx,dy,complex)
%
% Shift an image by the specified dx and dy offsets. These offsets may be
% non-integer as the shift is implemented by introducing a Fourier domain
% tilt term.
%
% Parameters
% ----------
% image : m x n array
%   Image to be shifted
%
% dx, dy : float or 1 x n arrays
%   (x, y) or (column, row) shift. If dx and dy are vectors, the result will
%   be a shift-and-add over the set of specified shifts. This is useful for
%   simulating jitter effects or motion blur.
%
% complex : bool, optional
%   If true, the complex result is returned. If false (default), only the
%   real portion of thre result is returned.
%
% Returns
% -------
% result : m x n array
%   Shifted image
%


if nargin < 4
    complex = false;
end

if (dx==0) && (dy==0)
    result = image;
    return;
end

if isempty(dx) || isempty(dy)
    result = image;
    return;
end

[m, n] = size(image);

N = length(dx);

xc = floor(n/2)+1;
yc = floor(m/2)+1;

[X, Y] = meshgrid(1:n,1:m);

X = X - xc;
Y = Y - yc;

%I = fftshift(fft2(image));
I = fft2(image);

T = 0*I;

for k=1:N
    
    px = -2*pi*(X)/n * dx(k);
    py = -2*pi*(Y)/m * dy(k);
    
    T = T + exp(1i * (px + py));    
end

%I = I .* T;
I = I .* fftshift(T);

%result = real(ifft2(fftshift(I)));
if complex
    result = ifft2(I);
else
    result = real(ifft2(I));
end

end
function vprintf(verbose, varargin)
% VPRINTF passes fprintf arguments if verbose is true.
%
% vprintf(verbose, <fprintf args>)
%
% Parameters
% ----------
% verbose : bool
%   If true, passes calls fprintf with arguments in varargin
%
% See Also
% --------
% FPRINTF
%


if verbose
    fprintf(varargin{:})
end
end

function [psf] = opd2psf(pupil,opd,lambda)

      P = pupil .* exp(i * 2*pi/lambda * opd);

      psf = fftshift(fft2(P));

      psf = psf .* conj(psf);

      psf = psf ./ sum(psf(:));

end
function [c] = circ(n,r)
    c= zeros(n);
    [x,y] = meshgrid(-n/2:n/2-1,-n/2:n/2-1);
    x = x+0.5;
    y = y+0.5;
    z = sqrt(x.^2+y.^2);
    c(z<r)=1;
end



                                                                                             