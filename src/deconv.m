function [x,bkg] = deconv(cfg,A,AAIinv,y,x,regularization)
%TODO: `eval(regularization)` is a very bad SWing choice!
%      better would be to create a class for ADMM, create methods
%      deconv/debias and deal with this internaly, if necessary
    fnAdmm = str2func(['ADMM_',cfg.noiseModel]);
    [x,bkg] = fnAdmm(cfg,A,AAIinv,y,x,regularization);
end

function [x,bkg] = ADMM_G(cfg,A,AAIinv,y,x,regularization)

    active_set = (x > 0);   % `active_set` is required for debias regularizer
    
    [bkg,w] = estimBkgAndWeights(A,GPU.from(y),cfg.bkg.w_factor,cfg.bkg.blur_sigma,cfg.bkg.estim_iters);   % `active_set` is required for deconv regularizer
    bkg = GPU.to(bkg);
    w = GPU.to(w);
    
    u2 = x; u3 = x;
    d2 = GPU.to(zeros(size(u2))); d3 = GPU.to(zeros(size(u3)));
    k = 0;
    while k < cfg.max_iter
        ksi2 = u2 + d2; ksi3 = u3 + d3;
        gamma = A.corr(y) + cfg.mu_regularizer*ksi2 + cfg.mu_positive*ksi3;
        x = AAIinv * gamma;
        nu2 = x - d2; nu3 = x - d3;
        u2 = eval(regularization);  % regularization is applied on `nu2`
        u3 = max(nu3,0);
        d2 = d2 - (x - u2); d3 = d3 - (x - u3);
        k = k + 1;
        
        % background update
        if (mod(k,cfg.bkg_update) == 0)
            [bkg,w] = estimBkgAndWeights(A,GPU.from(y - A.conv(u3)),cfg.bkg.w_factor,cfg.bkg.blur_sigma,cfg.bkg.update_iters);
            bkg = GPU.to(bkg);
            w = GPU.to(w);
        end
        
        % debug info
        if((cfg.debug > 0) && (mod(k,cfg.debug)==0))
            showProgress(k,cfg.z_range,A,y,u3,y,bkg);
        end
        
    end
    x = u3;
    if(cfg.debug > 0)
        est = A.conv(x) + bkg;
        fprintf('Final: l2err = %f\n',norm(y(:) - est(:)));
    end
    
end

function [x,bkg] = ADMM_P(cfg,A,AAIinv,y,x,regularization)

    active_set = (x > 0);   % `active_set` is required for debias regularizer
    
    [bkg,w] = estimBkgAndWeights(A,GPU.from(y),cfg.bkg.w_factor,cfg.bkg.blur_sigma,cfg.bkg.estim_iters);   % `active_set` is required for deconv regularizer
    bkg = GPU.to(bkg);
    w = GPU.to(w);
    
    u1 = y; u2 = x; u3 = x;
    d1 = GPU.to(zeros(size(u1))); d2 = GPU.to(zeros(size(u2))); d3 = GPU.to(zeros(size(u3)));
    k = 0;
    while k < cfg.max_iter
        ksi1 = u1 + d1; ksi2 = u2 + d2; ksi3 = u3 + d3;
        gamma = A.corr(ksi1) + cfg.mu_regularizer*ksi2 + cfg.mu_positive*ksi3;
        x = AAIinv * gamma;
        yhat = A.conv(x) + bkg;
        nu1 = yhat - d1; nu2 = x - d2; nu3 = x - d3;
        u1 = (nu1 - 1./cfg.mu + sqrt(max((nu1 - 1./cfg.mu).^2 + 4.*y./cfg.mu,0))) ./ 2;
        u2 = eval(regularization);  % regularization is applied on `nu2`
        u3 = max(nu3,0);
        d1 = d1 - (yhat - u1); d2 = d2 - (x - u2); d3 = d3 - (x - u3);
        k = k + 1;
        
        % background update
        if (mod(k,cfg.bkg_update) == 0)
            [bkg,w] = estimBkgAndWeights(A,GPU.from(y - A.conv(u3)),cfg.bkg.w_factor,cfg.bkg.blur_sigma,cfg.bkg.update_iters);
            bkg = GPU.to(bkg);
            w = GPU.to(w);
        end
        
        % debug info
        if((cfg.debug > 0) && (mod(k,cfg.debug)==0))
            showProgress(k,cfg.z_range,A,y,u3,u1,bkg);
        end
        
    end
    x = u3;
    if(cfg.debug > 0)
        est = A.conv(x) + bkg;
        fprintf('Final: l2err = %f\n',norm(y(:) - est(:)));
    end
    
end

function tx = soft(x,thr)
    tx = sign(x) .* max(abs(x) - thr, 0);
end

function fx = fixed_supp(x,supp)
    fx = GPU.to(zeros(size(x)));
    fx(supp) = x(supp);
end

function showProgress(iter,z_range,A,Y,Xest,Yest,Best)
    figure(1);
    zslices = length(z_range);
    for zi=1:zslices
        subplot(ceil(sqrt(zslices)),ceil(sqrt(zslices)),zi),imagesc(Xest(:,:,zi)),colorbar,title(sprintf('%f nm; iter=%d',z_range(zi),iter)),drawnow;
    end
    Yhat = A.conv(Xest) + Best;
    fprintf('l2err = %g\n',norm(Y(:) - Yhat(:)));
    figure(2);
    subplot(2,3,1),imagesc(imfirst(Y)),title('orig (y)'),colorbar,drawnow;
    subplot(2,3,2),imagesc(imfirst(Best)),title('bkg'),colorbar,drawnow;
    subplot(2,3,3),imagesc(imfirst(Yhat)),title('estim'),colorbar,drawnow;
    subplot(2,3,4),imagesc(imfirst(Y-Best-Yhat)),title('residual'),colorbar,drawnow;
    subplot(2,3,5),imagesc(imfirst(Y-Best)),title('orig-bkg'),colorbar,drawnow;
    subplot(2,3,6),imagesc(imfirst(Yest)),title('u1 (new y)'),colorbar,drawnow;
    if size(Y,3) == 2   % biplane
        figure(3);
        subplot(2,3,1),imagesc(imsecond(Y)),title('orig (y)'),colorbar,drawnow;
        subplot(2,3,2),imagesc(imsecond(Best)),title('bkg'),colorbar,drawnow;
        subplot(2,3,3),imagesc(imsecond(Yhat)),title('estim'),colorbar,drawnow;
        subplot(2,3,4),imagesc(imsecond(Y-Best-Yhat)),title('residual'),colorbar,drawnow;
        subplot(2,3,5),imagesc(imsecond(Y-Best)),title('orig-bkg'),colorbar,drawnow;
        subplot(2,3,6),imagesc(imsecond(Yest)),title('u1 (new y)'),colorbar,drawnow;
    end
end

function I = imfirst(im)
    I = im(:,:,1);
end

function I = imsecond(im)
    I = im(:,:,2);
end