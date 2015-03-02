function [cfg,file] = default_config(fpath,camera,calibration,mu_regularizer,reg_mult)

    %% Top-level settings
    file.path = fpath;
    [file.info,file.frames,file.width,file.height] = IO.getImageInfo(fpath,calibration);
    
    cfg.progressbar = true;
    cfg.calibration = calibration;
    cfg.camera = camera;
    
    cfg.thrIrel = 0.1;
    cfg.thrIabs = 10.0;
    cfg.noiseModel = 'P';   %'G' or 'P' for either Gaussian, or Poisson
    
    %% Background estimation
    cfg_bkg.w_factor = reg_mult / mu_regularizer;
    cfg_bkg.blur_sigma = 10.0;
    cfg_bkg.estim_iters = 10;
    cfg_bkg.update_iters = 1;
    
    %% ADMM (l1-regularized)
    cfg.deconv.debug = 0;
    cfg.deconv.max_iter = 500;
    cfg.deconv.bkg_update = 1;
    cfg.deconv.zoom = 3;
    cfg.deconv.z_step = 100;
    cfg.deconv.z_range = -400:cfg.deconv.z_step:+400;
    cfg.deconv.mu = 1e-3;
    cfg.deconv.mu_regularizer = mu_regularizer;
    cfg.deconv.mu_positive = 1.0;
    cfg.deconv.noiseModel = cfg.noiseModel;
    cfg.deconv.convEval = 'fft';%{'fft','zero','replicate','circular','symmetric'}
    cfg.deconv.hessMult = 'fft';%{'fft','zero','replicate','circular','symmetric'}
    cfg.deconv.bkg = cfg_bkg;
    
    %% ADMM (fixed support)
    cfg.debias.debug = 0;
    cfg.debias.max_iter = 200;
    cfg.debias.bkg_update = 1;
    cfg.debias.mu = cfg.deconv.mu;
    cfg.debias.mu_regularizer = mu_regularizer;
    cfg.debias.mu_positive = cfg.deconv.mu_positive;
    cfg.debias.noiseModel = cfg.noiseModel;
    cfg.debias.zoom = cfg.deconv.zoom;
    cfg.debias.z_step = cfg.deconv.z_step;
    cfg.debias.z_range = cfg.deconv.z_range;
    cfg.debias.bkg = cfg_bkg;
    
    %% LM-fitting
    cfg.refine.debug = 0;
    cfg.refine.max_iter = 200;
    cfg.refine.bkg_update = 1;
	cfg.refine.psfEval = 'Analytic';%{'Analytic','Taylor1','Taylor2','Spline'}
	cfg.refine.jacobianEval = 'Analytic';%{'Analytic','Taylor1','Spline'}
    cfg.refine.noiseModel = cfg.noiseModel;
    cfg.refine.zoom = 3;
    cfg.refine.z_step = 10;
    cfg.refine.z_range = -500:cfg.refine.z_step:+500;
    cfg.refine.calibration = cfg.calibration;
    cfg.refine.peaks_min_distance = 3;%floor(1.86 * A.zoom * mean([mean(A.psf_sigmas_x),mean(A.psf_sigmas_y)]) / 2);
    cfg.refine.fitradius = max(file.width,file.height);%10;%7;%round(3*max([A.psf_sigmas_x(:);A.psf_sigmas_y(:)]));
    cfg.refine.lambda = 0.001;
    cfg.refine.max_lambda = 1e10;
    cfg.refine.min_improvement = 1e-4;
    cfg.refine.min_mol_update = 1e-6;
    cfg.refine.max_improvement_full_step = 1e-3;
    cfg.refine.lambda_factor = 10;
    cfg.refine.bkg = cfg_bkg;

end