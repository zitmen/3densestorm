function calcJacobian(is_biplane)
% run: calcJacobian(true);
%      calcJacobian(false);
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
    
    %mols = [16,16,0,100,0];
    %{
    mols = [ 0, 0,   0,100,0; ...
            32,16,-400,100,0; ...
            16,16,+400,100,0];
    %}
    %
    mols = [ 0.15, 0.18, -52,100,0; ...
            32.16,16.32,-356,100,0; ...
            16.11,16.88,+286,100,0];
    %
    fprintf('\n!!! TODO: resolve `dz` !!! the following are not really error, but actual values!\n\n');
    
    fprintf('Times:\n');    
    fprintf('Analytic:\n'); aJ = imshape(cfg,testAnalytic(psf,cfg,xgrid,ygrid,mols));
    fprintf('Taylor1:\n'); t1J = imshape(cfg,testTaylor1(psf,cfg,xgrid,ygrid,mols));
    fprintf('Spline:\n'); sJ = imshape(cfg,testSpline(psf,cfg,xgrid,ygrid,mols));
    
    err_t1J = t1J;% - aJ;
    err_sJ = sJ;% - aJ;
    
    for ii=1:size(aJ,4)  % planes
        for vi=3%1:size(aJ,3)      % variables {x,y,z,i,o}
            figure(vi)
            subplot(size(aJ,4),3,1+3*(ii-1)),imagesc(aJ(:,:,vi,ii)),title('Analytic image'),colorbar,drawnow
            subplot(size(aJ,4),3,2+3*(ii-1)),imagesc(err_t1J(:,:,vi,ii)),title(sprintf('Taylor1 approx. SSE: %g', sum(sum(err_t1J(:,:,vi,ii).^2)))),colorbar,drawnow
            subplot(size(aJ,4),3,3+3*(ii-1)),imagesc(err_sJ(:,:,vi,ii)),title(sprintf('Spline approx. SSE: %g', sum(sum(err_sJ(:,:,vi,ii).^2)))),colorbar,drawnow
        end
    end
    
    fprintf('\nPrecision:\n');
    fprintf('Taylor1 approx. SSE: %g\n',sum(err_t1J(:).^2));
    fprintf('Spline approx. SSE: %g\n',sum(err_sJ(:).^2));
end

function J = testAnalytic(psf,cfg,xgrid,ygrid,mols)
    fprintf('  Init...');
    tic;
    A = PSFrefine(psf,cfg.width,cfg.height,cfg.zoom,cfg.z_range,cfg.calibration,'Analytic','Analytic',cfg.noise);
    toc
    
    fprintf('  Eval...');
    tic;
    for ii=1:size(mols,1)
        if ii == 1
            J = A.calcJacobian(xgrid,ygrid,mols(ii,:),transform(mols(ii,:)));
        else
            J = J + A.calcJacobian(xgrid,ygrid,mols(ii,:),transform(mols(ii,:)));
        end
    end
    toc
end

function J = testTaylor1(psf,cfg,xgrid,ygrid,mols)
    fprintf('  Init...');
    tic;
    A = PSFrefine(psf,cfg.width,cfg.height,cfg.zoom,cfg.z_range,cfg.calibration,'Analytic','Taylor1',cfg.noise);
    toc
    
    fprintf('  Eval...');
    tic;
    for ii=1:size(mols,1)
        if ii == 1
            J = A.calcJacobian(xgrid,ygrid,mols(ii,:),transform(mols(ii,:)));
        else
            J = J + A.calcJacobian(xgrid,ygrid,mols(ii,:),transform(mols(ii,:)));
        end
    end
    toc
end

function J = testSpline(psf,cfg,xgrid,ygrid,mols)
    fprintf('  Init...');
    tic;
    A = PSFrefine(psf,cfg.width,cfg.height,cfg.zoom,cfg.z_range,cfg.calibration,'Analytic','Spline',cfg.noise);
    toc
    
    fprintf('  Eval...');
    tic;
    for ii=1:size(mols,1)
        if ii == 1
            J = A.calcJacobian(xgrid,ygrid,mols(ii,:),transform(mols(ii,:)));
        else
            J = J + A.calcJacobian(xgrid,ygrid,mols(ii,:),transform(mols(ii,:)));
        end
    end
    toc
end

function I = imshape(cfg,J)
    planes = 1;
    if cfg.calibration.is_biplane, planes = 2; end;
    I = zeros(cfg.height,cfg.width,size(J,2),planes);
    imcnt = cfg.width*cfg.height;
    for ii=1:planes
        for vi=1:size(J,2)
            I(:,:,vi,ii) = reshape(J((ii-1)*imcnt+1:ii*imcnt,vi),[cfg.height,cfg.width]);
        end
    end
end

% copy-pasted from `refine.m`
function tmol = transform(mol)
    global ind nparams;
    tmol = mol;
    tmol(ind.I:nparams:end) = mol(ind.I:nparams:end).^2;
    tmol(ind.O:nparams:end) = mol(ind.O:nparams:end).^2;
end