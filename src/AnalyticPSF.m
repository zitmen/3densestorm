classdef AnalyticPSF
    
    methods (Static = true, Access = public)
        
        function psf = generate(low_w,low_h,zoom,z_range,cal)
            x_range = (1:1/zoom:low_w) - 0.5 - low_w/2;
            y_range = (1:1/zoom:low_h) - 0.5 - low_h/2;
            px = cal.px / zoom;
            center = 1/zoom/2;  % shifted a little to keep the peak in the middle of a pixel! if would be correct if `ceter=0`, BUT it would create a plateau of 4 neighbour pixels resulting in worse approximation!
            if cal.is_biplane
                psf = AnalyticPSF.genBiplane([x_range',x_range'],[y_range',y_range'],z_range,center,center,cal,px);
            else
                psf = AnalyticPSF.genAstigmatic(x_range,y_range,z_range,center,center,cal,px);
            end
        end
        
        function psf = genAstigmatic(x,y,z,x0,y0,cal,px)
        % x ..... x-range, e.g., -5:0.5:+5 [px]
        % y ..... y-range, e.g., -5:0.5:+5 [px]
        % z ..... z-range, e.g., -400:100:+400 [nm]
        % cal ... struct { w0 [nm], cx [nm], cy [nm], d [nm], fi [rad], px [nm] }
        %
        % returns struct { px [nm], im }
            [wx,wy] = AnalyticPSF.defocusGaussian(cal,z);
            z_step = 0; if length(z) > 1, z_step = abs(z(1)-z(2)); end;
            psf.voxel = [px,px,z_step];
            psf.stack = zeros(length(y),length(x),length(z));
            [xgrid,ygrid] = meshgrid(x,y);
            for zi = 1:length(z)
                psf.stack(:,:,zi) = AnalyticPSF.genGaussian(xgrid,ygrid,x0,y0,wx(zi),wy(zi),cal.fi);
            end
        end
        
        function psf = genBiplane(x,y,z,x0,y0,cal,px)
        % x ..... x-range, e.g., [[-5:0.5:+5]',[-5:0.5:+5]'] [px]
        % y ..... y-range, e.g., [[-5:0.5:+5]',[-5:0.5:+5]'] [px]
        % z ..... z-range, e.g., -400:100:+400 [nm]
        % cal ... struct { w0 [nm], {cx1,cx2} [nm], {cy1,cy2} [nm], d [nm], fi [rad], px [nm] }
        %
        % returns struct { px [nm], im }
            [wx,wy] = AnalyticPSF.defocusGaussian(cal,z);
            z_step = 0; if length(z) > 1, z_step = abs(z(1)-z(2)); end;
            psf.voxel = [px,px,z_step];
            psf.stack = zeros(size(y,1),size(x,1),length(z),2);
            planes = size(y,2);
            xgrid = zeros(size(y,1),size(x,1),planes);
            ygrid = zeros(size(y,1),size(x,1),planes);
            for ii = 1:planes
                [xgrid(:,:,ii),ygrid(:,:,ii)] = meshgrid(x(:,1),y(:,1));
            end
            for zi = 1:length(z)
                for ii = 1:planes
                    psf.stack(:,:,zi,ii) = AnalyticPSF.genGaussian(xgrid(:,:,ii),ygrid(:,:,ii),x0,y0,wx(ii,zi),wy(ii,zi),cal.fi) ./ planes;
                end
            end
        end
        
        %==================================================================================%
        function psf = genGaussian(xgrid,ygrid,x0,y0,xsigma,ysigma,fi)
            sinfi = sin(fi);
            cosfi = cos(fi);
            dx = ((xgrid - x0).*cosfi - (ygrid - y0).*sinfi);
            dy = ((xgrid - x0).*sinfi + (ygrid - y0).*cosfi);
            psf = exp(-0.5 .* ((dx./xsigma).^2 + (dy./ysigma).^2)) ./ (2*pi*xsigma*ysigma);
        end
        
        function [wx,wy] = defocusGaussian(cal,z)
            wx = zeros(length(cal.cx),length(z));
            wy = zeros(length(cal.cy),length(z));
            for zi=1:length(z)
                wx(:,zi) = (cal.w0./2) .* sqrt(1 + ((z(zi) - cal.cx)./cal.d).^2);
                wy(:,zi) = (cal.w0./2) .* sqrt(1 + ((z(zi) - cal.cy)./cal.d).^2);
            end
        end
        
    end
    
end

