function constructor(is_biplane)
% run: constructor(true);
%      constructor(false);
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
    
    cfg.convEval = 'fft';%{'fft','imfilter','call'}
    
    fprintf('matchingZ SSE: %g\n', matchingZ(cfg));
    fprintf('oversampledZ SSE: %g\n', oversampledZ(cfg));
    fprintf('undersampledZ SSE: %g\n', undersampledZ(cfg));
end

function SSE = matchingZ(cfg)
    psf = AnalyticPSF.generate(cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step:+500,cfg.calibration);
    psf.voxel = [80/cfg.zoom,80/cfg.zoom,cfg.z_step];
    
    A = PSFdeconv(psf,size(psf.stack,2),size(psf.stack,1),cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step:+500,cfg.calibration,cfg.convEval);
    
    for zz=1:size(psf.stack,3)
        psfNorm(:,:,zz,:) = psf.stack(:,:,zz,:) ./ sum(sum(sum(psf.stack(:,:,zz,:)))) .* cfg.zoom^2;
    end
    
    err = A.psf - psfNorm;
    SSE = sum(err(:).^2);
end

function SSE = oversampledZ(cfg)
    psf = AnalyticPSF.generate(cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step/2:+500,cfg.calibration);
    psf.voxel = [80/cfg.zoom,80/cfg.zoom,cfg.z_step/2];
    
    A = PSFdeconv(psf,size(psf.stack,2),size(psf.stack,1),cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step:+500,cfg.calibration,cfg.convEval);
    
    for zz=1:size(psf.stack,3)
        psfNorm(:,:,zz,:) = psf.stack(:,:,zz,:) ./ sum(sum(sum(psf.stack(:,:,zz,:)))) .* cfg.zoom^2;
    end
    
    err = A.psf - psfNorm(:,:,1:2:end,:);
    SSE = sum(err(:).^2);
end

function SSE = undersampledZ(cfg)
    psf = AnalyticPSF.generate(cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step:+500,cfg.calibration);
    psf.voxel = [80/cfg.zoom,80/cfg.zoom,cfg.z_step];
    
    A = PSFdeconv(psf,size(psf.stack,2),size(psf.stack,1),cfg.width,cfg.height,cfg.zoom,-500:cfg.z_step/2:+500,cfg.calibration,cfg.convEval);
    
    for zz=1:size(psf.stack,3)
        psfNorm(:,:,zz,:) = psf.stack(:,:,zz,:) ./ sum(sum(sum(psf.stack(:,:,zz,:)))) .* cfg.zoom^2;
    end
    
    err = A.psf(:,:,1:2:end,:) - psfNorm;
    SSE = sum(err(:).^2);
end