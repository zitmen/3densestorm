classdef PSFdeconv
% Conv/Corr eval:
%   fft....fast and preffered, but can consume a lot of memory, since
%          the size has to be identical to the size of zoomedIn analyzed
%          image; also if a measured PSF is used, if often has background
%          which causes problems and PSF needs to be padded and apodized
%   otherwise Matlab's `conv2` is called, but it's still possible to choose
%   boundary what conditions to use: zero,symmetric,circular,replicate
    
    properties (SetAccess = private)
        psf_bspline; psf; fpsf;
        x_range; y_range; z_range;
        hi_w; hi_h; hi_z;
        padd_w; padd_h;
        zoom;
        fnConv; fnCorr;
    end
    
    %% Internal padding of an image in case of non-fft operations with non-zero boundary conditions
    methods (Access = public)
        
        function P = paddImg(T,D,type)
            P = padarray(D,[T.padd_h,T.padd_w,0,0],type,'both');
        end
        
        function D = unpaddImg(T,P)
            D = P(1+T.padd_h:end-T.padd_h,1+T.padd_w:end-T.padd_w,:,:);
        end
        
        function P = paddZoomedImg(T,D,type)
            P = padarray(D,[T.padd_h*T.zoom,T.padd_w*T.zoom,0,0],type,'both');
        end
        
        function D = unpaddZoomedImg(T,P)
            D = P(T.zoom*T.padd_h+1:end-T.zoom*T.padd_h,T.zoom*T.padd_w+1:end-T.zoom*T.padd_w,:,:);
        end
        
    end
    
    %% Convolution: A*x
    methods (Access = protected)
        
        function M = conv_fft(T,U)
            M = GPU.to(zeros(size(U,1),size(U,2),size(T.psf,4)));
            for ii = 1:size(T.psf,4)
                M(:,:,ii) = T.unpadd(sum(real(ifft2( T.fpsf(:,:,:,ii) .* fft2(T.padd(U)) )),3));
            end
            M = T.zoomOut(M);
        end
        
        function M = conv_core(T,U)
            M = GPU.to(zeros(size(U,1),size(U,2),size(T.psf,4)));
            for ii = 1:size(T.psf,4)
                for zz = 1:size(T.psf,3)
                    M(:,:,ii) = M(:,:,ii) + conv2(U(:,:,zz),T.psf(:,:,zz,ii),'same');
                end
            end
            M = T.zoomOut(M);
        end
        
        function M = conv_zero(T,U)
            M = T.conv_core(U);
        end
        
        function M = conv_symmetric(T,U)
            M = T.unpaddImg(T.conv_core(T.paddZoomedImg(U,'symmetric')));
        end
        
        function M = conv_circular(T,U)
            M = T.unpaddImg(T.conv_core(T.paddZoomedImg(U,'circular')));
        end
        
        function M = conv_replicate(T,U)
            M = T.unpaddImg(T.conv_core(T.paddZoomedImg(U,'replicate')));
        end
        
    end
    
    %% Correlation: A'*y
    methods (Access = protected)
        
        function M = corr_fft(T,U)
            U = T.zoomIn(U);
            M = GPU.to(zeros(size(U,1),size(U,2),size(T.psf,3)));
            for zz = 1:size(T.psf,3)
                M(:,:,zz) = T.unpadd(sum(real(ifft2( conj(squeeze(T.fpsf(:,:,zz,:))) .* fft2(T.padd(U)) )),3));
            end
        end
        
        function M = corr_core(T,U)
            U = T.zoomIn(U);
            M = GPU.to(zeros(size(U,1),size(U,2),size(T.psf,3)));
            for zz = 1:size(T.psf,3)
                for ii = 1:size(T.psf,4)
                    M(:,:,zz) = M(:,:,zz) + conv2(U(:,:,ii),rot90(T.psf(:,:,zz,ii),2),'same');
                end
            end
        end
        
        function M = corr_zero(T,U)
            M = T.corr_core(U);
        end
        
        function M = corr_symmetric(T,U)
            M = T.unpaddZoomedImg(T.corr_core(T.paddImg(U,'symmetric')));
        end
        
        function M = corr_circular(T,U)
            M = T.unpaddZoomedImg(T.corr_core(T.paddImg(U,'circular')));
        end
        
        function M = corr_replicate(T,U)
            M = T.unpaddZoomedImg(T.corr_core(T.paddImg(U,'replicate')));
        end
        
    end
    
    %% Public methods
    methods (Access = public)
        
        function T = PSFdeconv(psf,psf_w,psf_h,low_w,low_h,zoom,z_range,cal,convEval)
        % Note: z-sampling of PSF is not based on z-coordinates but on sampling
        %       --> this can cause problems if the ranges aren't the same!!
        %           this could be fixed by adding an offset to the first value
            % initializing properties
            T.fnConv = str2func(['conv_',convEval]);
			T.fnCorr = str2func(['corr_',convEval]);
            T.zoom = zoom;
            T.hi_w = (low_w - 1) * T.zoom + 1;
            T.hi_h = (low_h - 1) * T.zoom + 1;
            T.hi_z = length(z_range);
            down_x = 1:T.zoom:T.hi_w;
            down_y = 1:T.zoom:T.hi_h;
            % B-spline approximation of PSF
            xy_step = cal.px/zoom;
            z_step = psf.voxel(3);
            if(length(z_range) > 1), z_step = abs(z_range(2) - z_range(1)); end;
            T.x_range = psf.voxel(2) + (xy_step .* ((1:psf_w)-1));
            T.y_range = psf.voxel(1) + (xy_step .* ((1:psf_h)-1));
            T.z_range = psf.voxel(3) + ( z_step .* ((1:length(z_range))-1));
            [xgrid,ygrid,zgrid] = meshgrid(T.x_range,T.y_range,T.z_range);
            for ii = 1:size(psf.stack,4)
                if size(psf.stack,3) > 1    % 3D
                    T.psf_bspline{ii} = bsarray(psf.stack(:,:,:,ii),'degree',3,'elementSpacing',psf.voxel,'lambda',0);
                    psf_stack(:,:,:,ii) = interp3(T.psf_bspline{ii},xgrid,ygrid,zgrid,0);
                else   % 2D
                    T.psf_bspline{ii} = bsarray(psf.stack(:,:,:,ii),'degree',3,'elementSpacing',psf.voxel(1:2),'lambda',0);
                    psf_stack(:,:,:,ii) = interp2(T.psf_bspline{ii},xgrid,ygrid,0);
                end
            end
            % PSF normalization
            T.psf = psf_stack;
            for zz=1:size(T.psf,3)
                T.psf(:,:,zz,:) = T.psf(:,:,zz,:) ./ sum(sum(sum(T.psf(:,:,zz,:)))) .* T.zoom^2;
            end
            T.padd_w = floor(ceil(size(T.psf,2)/T.zoom) / 2);
            T.padd_h = floor(ceil(size(T.psf,1)/T.zoom) / 2);
            % FFT
            if strcmp(convEval,'fft')
                % pre-computation of OTF (padded if necessary to ensure periodicity!)
                T.padd_w = (T.zoom - 1) - (T.hi_w - down_x(end));
                T.padd_h = (T.zoom - 1) - (T.hi_h - down_y(end));
                for zz=1:size(T.psf,3)
                    for ii=1:size(T.psf,4)
                        T.fpsf(:,:,zz,ii) = fft2(ifftshift(T.padd(T.psf(:,:,zz,ii))));
                    end
                end
                T.fpsf = GPU.to(T.fpsf);
            end
            T.psf = GPU.to(T.psf);
        end
        
        function AA2Iinv = getAAIinv(T,mu1,mu2,hessMult)
            % pre-computation of OTF (not padded!)
            if T.padd_w == 0 & T.padd_h == 0
                fpsf = T.fpsf;
            else
                for zz=1:size(T.psf,3)
                    for ii=1:size(T.psf,4)
                        fpsf(:,:,zz,ii) = fft2(ifftshift(T.psf(:,:,zz,ii)));
                    end
                end
                fpsf = GPU.to(fpsf);
            end
            % calculate A'*A
            AtA = cell(size(T.psf,3),size(T.psf,3));
            for z = 1:size(T.psf,3)
                for z2 = 1:size(T.psf,3)
                    AtA{z,z2} = zeros(size(fpsf,1),size(fpsf,2));
                    for ii = 1:size(T.psf,4)
                        AtA{z,z2} = AtA{z,z2} + conj(fpsf(:,:,z,ii)) .* fpsf(:,:,z2,ii);
                    end
                    AtA{z,z2} = AtA{z,z2} ./ T.zoom^2;
                    if z==z2
                        AtA{z,z2} = AtA{z,z2} + mu1 + mu2;
                    end
                end
            end
            % calculate inverse of the regularized Hessian
            AA2Iinv = PSF2Inv(T,AtA,hessMult);
        end
        
        function M = conv(T,U)
            M = T.fnConv(T,U);
        end
        
        function M = corr(T,U)
            M = T.fnCorr(T,U);
        end
        
        function Z = zoomIn(T,D)
            Z = GPU.to(zeros((size(D,1)-1)*T.zoom+1,(size(D,2)-1)*T.zoom+1,size(D,3)));
            Z(1:T.zoom:end,1:T.zoom:end,:) = D;
        end
        
        function D = zoomOut(T,Z)
            D = Z(1:T.zoom:end,1:T.zoom:end,:);
        end
        
        function P = padd(T,D)
            P = padarray(D,[T.padd_h,T.padd_w,0,0],0,'post');
        end
        
        function D = unpadd(T,P)
            D = P(1:end-T.padd_h,1:end-T.padd_w,:,:);
        end
         
    end
    
end
