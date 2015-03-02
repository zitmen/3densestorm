addpath(genpath('../../src'));

run_in_parallel = true;

camera.pixelsize   = 100;    % nm
camera.offset      = 100;    % digital units
camera.photons2ADU = 2.0;    % ADU conversion rate
camera.gain        = 1.0;    % EM gain (multiplier)

calibration.is_biplane = false;
calibration.px =  camera.pixelsize;        % nm
calibration.w0 =  323.77 / calibration.px; % px
calibration.cx =  150.00;   % nm
calibration.cy =  150.00;   % nm
calibration.d  =  400.00;   % nm
calibration.fi =    0.00;   % rad

fname = 'tubulins';
[cfg,file] = default_config([fname,'.tif'],camera,calibration,0.1,15);
cfg.deconv.z_range = 0;
cfg.debias.z_range = 0;
cfg.refine.z_range = 0;
cfg.deconv.z_step = 0;
cfg.debias.z_step = 0;
cfg.refine.z_step = 0;
%file.frames = 1:100;

psf_deconv = AnalyticPSF.generate(file.width,file.height,cfg.deconv.zoom,cfg.deconv.z_range,cfg.calibration);
psf_deconv.voxel = [cfg.calibration.px/cfg.deconv.zoom,cfg.calibration.px/cfg.deconv.zoom,cfg.deconv.z_step];

psf_refine = AnalyticPSF.generate(file.width,file.height,cfg.refine.zoom,cfg.refine.z_range,cfg.calibration);
psf_refine.voxel = [cfg.calibration.px/cfg.refine.zoom,cfg.calibration.px/cfg.refine.zoom,cfg.refine.z_step];

denseSTORM(cfg,file,[fname,'-results.csv'],psf_deconv,psf_refine,run_in_parallel);