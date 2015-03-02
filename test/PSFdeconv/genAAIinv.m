function genAAIinv(is_biplane)
% run: genAAIinv(true);
%      genAAIinv(false);
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
    
    cfg.mu1 = 0.05;
    cfg.mu2 = 1;
    
    cfg.convEval = 'fft';
    cfg.corrEval = 'fft';
    cfg.hessMult = 'fft';
    
    psf = AnalyticPSF.generate(cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step:+500,cfg.calibration);
    psf.voxel = [80/cfg.zoom,80/cfg.zoom,cfg.z_step];
    A = PSFdeconv(psf,size(psf.stack,2),size(psf.stack,1),cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step:+500,cfg.calibration,cfg.convEval);
    
    fAi = A.getAAIinv(cfg.mu1,cfg.mu2,'fft');
    cAi = A.getAAIinv(cfg.mu1,cfg.mu2,'circular');
    
    y = A.zoomOut(psf.stack(:,:,5,:));
    Aty = A.corr(y);
    
    fprintf('Times:\n');
    fprintf('  fft.....'); tic; fx = fAi*Aty; toc
    fprintf('  conv2....'); tic; cx = cAi*Aty; toc
    
    fprintf('\nPrecision:\n');
    fprintf('  SSE: fft-conv2: %g\n',sum((fx(:)-cx(:)).^2));
end