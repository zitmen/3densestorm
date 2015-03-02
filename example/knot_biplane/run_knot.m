addpath(genpath('../../src'));

run_in_parallel = true;

camera.pixelsize   =  80;    % nm
camera.offset      = 100;    % digital units
camera.photons2ADU = 1.0;    % ADU conversion rate
camera.gain        = 1.0;    % EM gain (multiplier)

calibration.is_biplane = true;
calibration.divide_dim = 2;
calibration.px =  camera.pixelsize;        % nm
calibration.w0 =  218.40 / calibration.px; % px
calibration.cx =[-150;+150];   % nm
calibration.cy =[-150;+150];   % nm
calibration.d  =  400.00;      % nm
calibration.fi =    0.00;      % rad

fname = 'knot';
[cfg,file] = default_config([fname,'.tif'],camera,calibration,0.1,10);

psf_deconv = AnalyticPSF.generate(file.width,file.height,cfg.deconv.zoom,cfg.deconv.z_range,cfg.calibration);
psf_deconv.voxel = [cfg.calibration.px/cfg.deconv.zoom,cfg.calibration.px/cfg.deconv.zoom,cfg.deconv.z_step];

psf_refine = AnalyticPSF.generate(file.width,file.height,cfg.refine.zoom,cfg.refine.z_range,cfg.calibration);
psf_refine.voxel = [cfg.calibration.px/cfg.refine.zoom,cfg.calibration.px/cfg.refine.zoom,cfg.refine.z_step];

denseSTORM(cfg,file,[fname,'-results.csv'],psf_deconv,psf_refine,run_in_parallel);