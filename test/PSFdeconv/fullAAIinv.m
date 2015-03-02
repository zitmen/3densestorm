function fullAAIinv()
% run: fullAAIinv();
%      --> tests astigmatism only as the extension is straight forward
    addpath(genpath('../../src'));
    
    cfg.width = 15;
    cfg.height = 15;
    cfg.zoom = 3;
    cfg.z_step = 200;
    cfg.z_range = -200:cfg.z_step:+200;
    
    cfg.calibration.is_biplane = false;
    cfg.calibration.px =   80.00; % nm
    cfg.calibration.w0 =    1.50; % px
    cfg.calibration.d  =  400.00; % nm
    cfg.calibration.fi =    0.00; % rad
    cfg.calibration.cx =  150.00;  % nm
    cfg.calibration.cy = -150.00;  % nm
    
    cfg.mu1 = 0.1;
    cfg.mu2 = 1;
    
    method = 'fft';%{'fft','circular'}
    shift = 1; if strcmp(method,'fft'), shift = 0; end;
    
    cfg.convEval = method;
    cfg.corrEval = method;
    cfg.hessMult = method;
    
    intensity = 100;
    
    psf = AnalyticPSF.generate(cfg.width,cfg.height,cfg.zoom,cfg.z_range,cfg.calibration);
    psf.voxel = [80/cfg.zoom,80/cfg.zoom,cfg.z_step];
    A = PSFdeconv(psf,size(psf.stack,2),size(psf.stack,1),cfg.width,cfg.height,cfg.zoom,cfg.z_range,cfg.calibration,cfg.convEval);
    y = A.zoomOut(intensity .* psf.stack(:,:,1,:));
    
    fAAIinv = A.getAAIinv(cfg.mu1,cfg.mu2,cfg.hessMult);
    aAAIinv = genFullAAIinv(cfg,A);
    
    Aty = A.corr(y);
    fx = fAAIinv*Aty;
    ax = circshift(flip3(reshape(aAAIinv*vect(flip3(Aty)),size(fx))),[shift,shift]);
    
    subplot(2,3,1),imagesc(ax(:,:,1)),title('alg'),colorbar
    subplot(2,3,2),imagesc(ax(:,:,2)),title('alg'),colorbar
    subplot(2,3,3),imagesc(ax(:,:,3)),title('alg'),colorbar
    subplot(2,3,4),imagesc(fx(:,:,1)),title('fft'),colorbar
    subplot(2,3,5),imagesc(fx(:,:,2)),title('fft'),colorbar
    subplot(2,3,6),imagesc(fx(:,:,3)),title('fft'),colorbar
    
    fprintf('SSE(fft-alg): %g\n',sum((fx(:) - ax(:)).^2));
end

function AAIinv = genFullAAIinv(cfg,psf)

    A = genFullA(cfg,psf);
    AtA = A'*A;
    I = eye(size(AtA));
    AAI = AtA + cfg.mu1.*I + cfg.mu2.*I;
    AAIinv = inv(AAI);

end

function A3D = genFullA(cfg,psf)

    h = cell(length(cfg.z_range),1);
    for zi = 1:length(h)
        h{zi} = psf.psf(:,:,zi);
    end

    I = GPU.to(zeros(psf.hi_h,psf.hi_w));
    A = GPU.to(zeros(cfg.height*cfg.width, psf.hi_h*psf.hi_w));

    Az = cell(length(h),1);
    for zi=1:length(h)
        ii = 0;
        for c=1:size(I,2)
            for r=1:size(I,1)
                ii = ii + 1;
                I(r,c) = 1;
                
                Ih = conv2(I,h{zi},'same');
                tmp = psf.zoomOut(Ih);
                A(:,ii) = tmp(:);

                I(r,c) = 0;
            end
        end
        Az{zi} = A;
    end
    
    A3D = GPU.to(zeros(size(A,1), length(h)*size(A,2)));
    for zi=1:length(h)
        A3D(:,(zi-1)*size(A,2)+1:zi*size(A,2)) = Az{zi};
    end
    
end

function v = vect(M)
    v = M(:);
end

function F = flip3(M)
    F = M(:,:,end:-1:1);
end