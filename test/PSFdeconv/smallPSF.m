function smallPSF()
% run: smallPSF();
%      --> tests astigmatism only as the extension is straight forward
    addpath(genpath('../../src'));
    
    cfg.width = 32;
    cfg.height = 32;
    cfg.zoom = 3;
    cfg.z_step = 200;
    cfg.z_range = -200:cfg.z_step:+200;
    
    cfg.calibration.is_biplane = false;
    cfg.calibration.px =   80.00; % nm
    cfg.calibration.w0 =    2.73; % px
    cfg.calibration.d  =  400.00; % nm
    cfg.calibration.fi =    0.00; % rad
    cfg.calibration.cx =  150.00; % nm
    cfg.calibration.cy = -150.00; % nm
    
    cfg.method = 'zero';
    cfg.mu1 = 0.1;
    cfg.mu2 = 1;
    
    psf = AnalyticPSF.generate(cfg.width,cfg.height,cfg.zoom,cfg.z_range,cfg.calibration);
    psf.voxel = [80/cfg.zoom,80/cfg.zoom,cfg.z_step];
    
    zoom = cfg.zoom;%2;4;
    [ly0,ly1] = largeCall(cfg,zoom,psf);
    [sy0,sy1] = smallCall(cfg,zoom);
    
    fprintf('iter=0: SSE(large-small): %g\n',sum((ly0(:) - sy0(:)).^2));
    fprintf('iter=1: SSE(large-small): %g\n',sum((ly1(:) - sy1(:)).^2));
end

function [y,y1] = largeCall(cfg,zoom,psf)
% size of PSF == size of image; fft/call/imfilter can be used in this case
    xy = cfg.width;
    A = PSFdeconv(psf,(xy-1)*zoom+1,(xy-1)*zoom+1,cfg.width,cfg.height,zoom,cfg.z_range,cfg.calibration,cfg.method);
    AAIinv = A.getAAIinv(cfg.mu1,cfg.mu2,cfg.method);
    
    x = getData(cfg.width,cfg.height,length(cfg.z_range),zoom);
    y = A.conv(x);
    subplot(2,2,1),imagesc(y),title(sprintf('large %s, iter=0',cfg.method)),colorbar
    
    Aty = A.corr(y);
    x1 = AAIinv*Aty;
    y1 = A.conv(x1);
    subplot(2,2,2),imagesc(y1),title(sprintf('large %s, iter=1',cfg.method)),colorbar
end

function [y,y1] = smallCall(cfg,zoom)
% size of PSF < size of image; call/imfilter can be used in this case (also applicable to PSF > size)
    psf = AnalyticPSF.generate(16,16,cfg.zoom,cfg.z_range,cfg.calibration);
    psf.voxel = [80/cfg.zoom,80/cfg.zoom,cfg.z_step];
    
    xy = 16;
    A = PSFdeconv(psf,(xy-1)*zoom+1,(xy-1)*zoom+1,cfg.width,cfg.height,zoom,cfg.z_range,cfg.calibration,cfg.method);
    AAIinv = A.getAAIinv(cfg.mu1,cfg.mu2,cfg.method);
    
    x = getData(cfg.width,cfg.width,length(cfg.z_range),zoom);
    y = A.conv(x);
    subplot(2,2,3),imagesc(y),title(sprintf('small %s, iter=0',cfg.method)),colorbar
    
    Aty = A.corr(y);
    x1 = AAIinv*Aty;
    y1 = A.conv(x1);
    subplot(2,2,4),imagesc(y1),title(sprintf('small %s, iter=1',cfg.method)),colorbar
end

function x = getData(w,h,z,zoom)
    x = GPU.to(zeros(h*zoom,w*zoom,z));
    x( 8*zoom,12*zoom,1) = 100;
    x(13*zoom, 8*zoom,2) = 100;
    x(26*zoom,25*zoom,3) = 100;
end