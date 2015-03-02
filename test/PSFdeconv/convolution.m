function convolution(is_biplane)
% run: convolution(true);
%      convolution(false);
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
    
    x = GPU.to(zeros(94,94,11));
    x(25,25, 1) = 100;
    x(45,45, 6) = 100;
    x(45,70,11) = 100;
    
    fI = fftConv(cfg,psf,x);
    cI = callConv(cfg,psf,x);
    
    for ii=1:size(fI,3)
        subplot(size(fI,3),3,(ii-1)*3+1),imagesc(fI(:,:,ii)),title('fft'),colorbar
        subplot(size(fI,3),3,(ii-1)*3+3),imagesc(cI(:,:,ii)),title('conv2'),colorbar
    end
    
    fprintf('SSE(fft-conv2): %g\n',sum((fI(:) - cI(:)).^2));
end

function I = fftConv(cfg,psf,x)
    A = PSFdeconv(psf,size(psf.stack,2),size(psf.stack,1),cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step:+500,cfg.calibration,'fft');
    I = A.conv(x);  % there is a slight shift (caused by padding which ensures correct periodicity for fft)
end

function I = callConv(cfg,psf,x)
    A = PSFdeconv(psf,size(psf.stack,2),size(psf.stack,1),cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step:+500,cfg.calibration,'circular');
    I = A.conv(x);
end