function denseSTORM(cfg,file,resfile,psf_deconv,psf_refine,parallel)
% 3denseSTORM
tic
    warning('off','all');
    fprintf('initialization...\n');
    
    if parallel && isempty(gcp('nocreate'))
        parpool('local');
        pctRunOnAll warning('off','all');
    end
    
    % init deconv admm
    A = PSFdeconv(psf_deconv,size(psf_deconv.stack,2),size(psf_deconv.stack,1),file.width,file.height,cfg.deconv.zoom,cfg.deconv.z_range,cfg.calibration,cfg.deconv.convEval);
    AAIinv = A.getAAIinv(cfg.deconv.mu_regularizer,cfg.deconv.mu_positive,cfg.deconv.hessMult);
    
    % init refinement
    refA = PSFrefine(psf_refine,size(psf_refine.stack,2),size(psf_refine.stack,1),cfg.refine.zoom,cfg.refine.z_range,cfg.calibration,cfg.refine.psfEval,cfg.refine.jacobianEval,cfg.refine.noiseModel);
toc
tic
    if cfg.progressbar, progress = ProgressBar(length(file.frames),'verbose',0); end;
    fprintf('Analyzing dataset consisting of %d frames\n',length(file.frames));
    header = {'"frame"','"x [nm]"','"y [nm]"','"z [nm]"','"I [photon]"','"offset [photon]"'};
    IO.writeHeader(resfile,header);
    startTime = now();
    if parallel
        parfor frame=file.frames
            fprintf('\nprocessing frame %d...\n',frame);
            analyzeImage(cfg,A,AAIinv,refA,file,frame,resfile);
            if cfg.progressbar, printProgress(progress.progress(),startTime); end;
        end
    else
        for frame=file.frames
            fprintf('\nprocessing frame %d...\n',frame);
            analyzeImage(cfg,A,AAIinv,refA,file,frame,resfile);
            if cfg.progressbar, printProgress(progress.progress(),startTime); end;
        end
    end
    fprintf('\n');
    if cfg.progressbar, progress.stop(); end;
toc
tic
    fprintf('sortng the results...\n');
    IO.sortResults(resfile,header);
    
    fprintf('\nstopping...\n');
    if parallel && ~isempty(gcp('nocreate'))
        delete(gcp('nocreate'));
    end
toc
end

function analyzeImage(cfg,A,AAIinv,refA,file,frame,resfile)
    % load a frame of sequence
    im = GPU.to(IO.readImage(file,cfg.calibration,frame));
    y = (im - cfg.camera.offset) .* cfg.camera.photons2ADU ./ cfg.camera.gain;    % convert to photons
    
    % deconvolution
    x = GPU.to(zeros(A.hi_h,A.hi_w,A.hi_z));
    x = deconv(cfg.deconv,A,AAIinv,y,x,'soft(nu2,w)');
    x(x < max(max(x(:))*cfg.thrIrel,cfg.thrIabs)) = 0;   % rel% of maximum intensity && thr >= abs
    if (sum(x(:) > 0) > 0)
        active_set = (x > 0);
        x = deconv(cfg.debias,A,AAIinv,y,x,'fixed_supp(nu2,active_set)');
        x(~active_set(:)) = 0;
    end
    thr = max(max(x(:))*cfg.thrIrel,cfg.thrIabs);   % rel% of maximum intensity && thr >= abs
    x(x < 0) = 0;
    
    % refinement
    global ind nparams clamp;
    nparams = 5;
    ind = struct('X',1,'Y',2,'Z',3,'I',4,'O',5);
    clamp = [1.0,1.0,100.0,500.0,5.0];
    
    mol = detect_peaks(GPU.from(x),thr,cfg.deconv.zoom,cfg.deconv.z_range);
    res = refine(cfg.refine,refA,y,mol);
    
    % store results
    if(~isempty(res))
        results = [frame.*ones(size(res.X(:))),res.X(:),res.Y(:),res.Z(:),res.I(:),res.O(:)];
        IO.appendResults(resfile,results);
    end
end

function printProgress(pdone,startTime)
    elapsedTime = now() - startTime;
    fprintf('done %3.2f%%, elapsed: %s, remaining: %s\n',pdone,...
        formatTime(elapsedTime),formatTime(elapsedTime / pdone * (100 - pdone)));
end

function str = formatTime(t)
    oneSec = datenum('00/01/0000 00:00:1','dd/mm/yyyy HH:MM:SS');
    tSeconds = floor(t / oneSec);
    tMinutes = floor(tSeconds / 60);
    tSeconds = tSeconds - tMinutes*60;
    tHours = floor(tMinutes / 60);
    tMinutes = tMinutes - tHours*60;
    if tHours > 0
        str = sprintf('%d hrs %d min %d sec',tHours,tMinutes,tSeconds);
    elseif tMinutes > 0
        str = sprintf('%d min %d sec',tMinutes,tSeconds);
    else
        str = sprintf('%d sec',tSeconds);
    end
end