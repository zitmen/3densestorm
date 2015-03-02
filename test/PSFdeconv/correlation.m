function correlation(is_biplane)
% run: correlation(true);
%      correlation(false);
    addpath(genpath('../../src'));
    
    cfg.width = 32;
    cfg.height = 32;
    cfg.zoom = 3;
    cfg.z_step = 100;
    
    cfg.calibration.is_biplane = is_biplane;
    cfg.calibration.divide_dim = 2;
    cfg.calibration.px =   80.00; % nm
    cfg.calibration.w0 =    2.73; % px
    cfg.calibration.d  =  400.00; % nm
    cfg.calibration.fi =    0.00; % rad
    if is_biplane
        cfg.calibration.cx = [-150;+150];   % nm
        cfg.calibration.cy = [-150;+150];   % nm
    else
        cfg.calibration.cx =  150.00;  % nm
        cfg.calibration.cy = -150.00;  % nm
    end
    
    psf = AnalyticPSF.generate(cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step:+500,cfg.calibration);
    psf.voxel = [80/cfg.zoom,80/cfg.zoom,cfg.z_step];
    
    img = AnalyticPSF.generate(cfg.width,cfg.height,1,+400,cfg.calibration);
    y = squeeze(img.stack);
    
    fI = fftCorr(cfg,psf,y);
    cI = callCorr(cfg,psf,y);
    
    fprintf('SSE(fft-conv2): %g\n',sum((fI(:) - cI(:)).^2));
end

function I = fftCorr(cfg,psf,y)
    A = PSFdeconv(psf,size(psf.stack,2),size(psf.stack,1),cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step:+500,cfg.calibration,'fft');
    I = circshift(A.corr(y),[-2,-2,0]);  % there is a slight shift (caused by padding which ensures correct periodicity for fft)
end

function I = callCorr(cfg,psf,y)
    A = PSFdeconv(psf,size(psf.stack,2),size(psf.stack,1),cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step:+500,cfg.calibration,'circular');
    I = A.corr(y);
end