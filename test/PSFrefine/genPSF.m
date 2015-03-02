function genPSF(is_biplane)
% run: genPSF(true);
%      genPSF(false);
    addpath(genpath('../../src'));
    
    cfg.width = 32;
    cfg.height = 32;
    cfg.zoom = 3;
    cfg.noise = 'P';
    
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
    
    cfg.psfEval = 'Analytic';%{'Analytic','Taylor1','Taylor2','Spline'}
    cfg.jacobianEval = 'Analytic';%{'Analytic','Taylor1','Spline'}
    z_step = 10;
    cfg.z_range = -500:z_step:+500;
    
    global ind nparams;
    nparams = 5;
    ind = struct('X',1,'Y',2,'Z',3,'I',4,'O',5);
    
    psf = AnalyticPSF.generate(cfg.width,cfg.height,cfg.zoom,cfg.z_range,cfg.calibration);
    psf.voxel = [80/cfg.zoom,80/cfg.zoom,z_step];
    
    [xgrid,ygrid] = meshgrid((1:cfg.width)-0.5,(1:cfg.height)-0.5,1:length(cfg.calibration.cx));
    
    mols = [ 0.15, 0.18, -52,100,0; ...
            32.16,16.32,-356,100,0; ...
            16.11,16.88,+286,100,0];
    
    fprintf('Times:\n');
    fprintf('Analytic:\n'); aG = testAnalytic(psf,cfg,xgrid,ygrid,mols);
    fprintf('Taylor1:\n'); t1G = testTaylor1(psf,cfg,xgrid,ygrid,mols);
    fprintf('Taylor2:\n'); t2G = testTaylor2(psf,cfg,xgrid,ygrid,mols);
    fprintf('Spline:\n'); sG = testSpline(psf,cfg,xgrid,ygrid,mols);
    
    err_t1G = t1G - aG;
    err_t2G = t2G - aG;
    err_sG = sG - aG;
    
    for ii=1:size(xgrid,3)
        figure(ii)
        subplot(2,2,1),imagesc(aG(:,:,ii)),title('Analytic image'),colorbar,drawnow
        subplot(2,2,2),imagesc(err_t1G(:,:,ii)),title(sprintf('Taylor1 approx. SSE: %g', sum(sum(err_t1G(:,:,ii).^2)))),colorbar,drawnow
        subplot(2,2,3),imagesc(err_t2G(:,:,ii)),title(sprintf('Taylor2 approx. SSE: %g', sum(sum(err_t2G(:,:,ii).^2)))),colorbar,drawnow
        subplot(2,2,4),imagesc(err_sG(:,:,ii)),title(sprintf('Spline approx. SSE: %g', sum(sum(err_sG(:,:,ii).^2)))),colorbar,drawnow
    end
    
    fprintf('\nPrecision:\n');
    fprintf('Taylor1 approx. SSE: %g\n',sum(err_t1G(:).^2));
    fprintf('Taylor2 approx. SSE: %g\n',sum(err_t2G(:).^2));
    fprintf('Spline approx. SSE: %g\n',sum(err_sG(:).^2));
end

function G = testAnalytic(psf,cfg,xgrid,ygrid,mols)
    fprintf('  Init...');
    tic;
    A = PSFrefine(psf,cfg.width,cfg.height,cfg.zoom,cfg.z_range,cfg.calibration,'Analytic','Analytic',cfg.noise);
    toc
    
    fprintf('  Eval...');
    tic;
    G = zeros(size(xgrid));
    for ii=1:size(mols,1)
        G = G + A.genPSF(xgrid,ygrid,mols(ii,:));
    end
    toc
end

function G = testTaylor1(psf,cfg,xgrid,ygrid,mols)
    fprintf('  Init...');
    tic;
    A = PSFrefine(psf,cfg.width,cfg.height,cfg.zoom,cfg.z_range,cfg.calibration,'Taylor1','Analytic',cfg.noise);
    toc
    
    fprintf('  Eval...');
    tic;
    G = zeros(size(xgrid));
    for ii=1:size(mols,1)
        G = G + A.genPSF(xgrid,ygrid,mols(ii,:));
    end
    toc
end

function G = testTaylor2(psf,cfg,xgrid,ygrid,mols)
    fprintf('  Init...');
    tic;
    A = PSFrefine(psf,cfg.width,cfg.height,cfg.zoom,cfg.z_range,cfg.calibration,'Taylor2','Analytic',cfg.noise);
    toc
    
    fprintf('  Eval...');
    tic;
    G = zeros(size(xgrid));
    for ii=1:size(mols,1)
        G = G + A.genPSF(xgrid,ygrid,mols(ii,:));
    end
    toc
end

function G = testSpline(psf,cfg,xgrid,ygrid,mols)
    fprintf('  Init...');
    tic;
    A = PSFrefine(psf,cfg.width,cfg.height,cfg.zoom,cfg.z_range,cfg.calibration,'Spline','Analytic',cfg.noise);
    toc
    
    fprintf('  Eval...');
    tic;
    G = zeros(size(xgrid));
    for ii=1:size(mols,1)
        G = G + A.genPSF(xgrid,ygrid,mols(ii,:));
    end
    toc
end
