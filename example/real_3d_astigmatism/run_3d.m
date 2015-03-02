addpath(genpath('../../src'));

run_in_parallel = true;

camera.pixelsize   =  65;    % nm
camera.offset      =  99;    % digital units
camera.photons2ADU = 0.6;    % ADU conversion rate
camera.gain        = 1.0;    % EM gain (multiplier)

calibration.is_biplane = false;
calibration.px = camera.pixelsize;

fname = 'tubulins';
[cfg,file] = default_config([fname,'.tif'],camera,calibration,0.1,10);
cfg.thrIrel = 0.05;
cfg.refine.psfEval = 'Taylor2';
cfg.refine.jacobianEval = 'Taylor1';
file.frames = 1:100;    % process subsequence

% reading a measured PSF
fcal.path = 'calibration.tif';
[fcal.info,fcal.frames,fcal.width,fcal.height] = IO.getImageInfo(fcal.path,calibration);

% preparing the PSF for deconvolution (needs apodization near boundaries due to the background and its frequency response)
cfg.deconv.z_step = 100;
cfg.deconv.z_range = -400:cfg.deconv.z_step:+400;
psf_deconv.voxel = [cfg.calibration.px/cfg.deconv.zoom,cfg.calibration.px/cfg.deconv.zoom,100];
psf_deconv.stack = zeros(1+(file.width-1)*cfg.deconv.zoom,1+(file.height-1)*cfg.deconv.zoom,length(cfg.deconv.z_range));
wr = tukeywin(fcal.height*cfg.deconv.zoom,0.5);
wc = tukeywin(fcal.width*cfg.deconv.zoom,0.5);
[maskr,maskc] = meshgrid(wr,wc);
W = maskr .* maskc;
padw = ceil((size(psf_deconv.stack,2) - fcal.width*cfg.deconv.zoom) / 2);
padh = ceil((size(psf_deconv.stack,1) - fcal.height*cfg.deconv.zoom) / 2);
minval = inf; for frame=1:10:81, minval = min(min(min(IO.readImage(fcal,calibration,frame))),minval); end;
zi = 0;
for frame=1:10:81
    zi = zi + 1;
    tmp = padarray(W .* imresize(double(IO.readImage(fcal,calibration,frame)-minval),cfg.deconv.zoom,'Method','bicubic'),[padh,padw,0],0,'both');
    psf_deconv.stack(:,:,zi) = tmp(1:size(psf_deconv.stack,1),1:size(psf_deconv.stack,2));
end

% preparing PSF for refinement
cfg.refine.zoom = 1;
cfg.refine.z_step = 10;
cfg.refine.z_range = -400:cfg.refine.z_step:+400;
psf_refine.voxel = [cfg.calibration.px,cfg.calibration.px,cfg.refine.z_step];
minval = inf; for frame=1:81, minval = min(min(min(IO.readImage(fcal,calibration,frame))),minval); end;
for zi=fcal.frames
    psf_refine.stack(:,:,zi) = double(IO.readImage(fcal,calibration,frame)-minval);
end
cfg.refine.fitregion = min(size(psf_refine.stack,1),size(psf_refine.stack,2));

% note: the attached PSF is a raw image of a fluorescent bead and it is
%       not a good PSF by any means; this is just an example; due to
%       its high background and noise, it's not used for refinement
cfg.refine.max_iter = 0;

denseSTORM(cfg,file,[fname,'-results.csv'],psf_deconv,psf_refine,run_in_parallel);
