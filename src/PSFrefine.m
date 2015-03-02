classdef PSFrefine
    
    properties (SetAccess = private)
        cal;
        psf_bspline;
        psf_dx_bspline; psf_dy_bspline; psf_dz_bspline;
        psf_dI; psf_dx; psf_dy; psf_dz;
        psf_dxx; psf_dxy; psf_dxz;
        psf_dyx; psf_dyy; psf_dyz;
        psf_dzx; psf_dzy; psf_dzz;
        low_w; low_h;
        hi_w; hi_h; hi_z;
        zoom;
        z_range; z_step;
        xy_range; xy_step;
        fnGenPsf; fnCalcJacobian; fnCalcDataDiffs;
    end
    
    %% Helper functions (calculation of indices and xy/z shifts)
    methods (Access = private)
        
        function [rows,cols,xd,yd] = calcIndicesXY(T,xgrid,ygrid,plane,mol)
            global ind;
            
            xd = xgrid(:,:,plane) + (T.low_w/2 - mol(ind.X));
            yd = ygrid(:,:,plane) + (T.low_h/2 - mol(ind.Y));
            cols = (xd(1,:) >= max(0.5-1/T.zoom,0) & xd(1,:) <= T.low_w-0.5);%xy_range +- zoom/2
            rows = (yd(:,1) >= max(0.5-1/T.zoom,0) & yd(:,1) <= T.low_h-0.5);%xy_range +- zoom/2
            
            xd = xd(rows,cols);
            yd = yd(rows,cols);
        end
        
        function [zi,zd] = calcIndexDiffZ(T,mol)
            global ind;
            
            [~,zi] = min((T.z_range - mol(ind.Z)).^2);
            zd = (mol(ind.Z) - T.z_range(zi)) / T.z_step;
            if (abs(zd) > 1.0), zi = -1; end;    % `z` out of range
        end
        
        function [xd,yd,rr,cc] = calcDiffs(T,rows,cols,x0,y0)
            % this is done so stupid, because of the sampling, e.g. 1/3
            tmp = (T.xy_range - x0).^2; xv = min(tmp); xi = find((tmp-xv)<1e-6,1,'first');
            xd = (x0 - T.xy_range(xi)) / T.xy_step;
            tmp = (T.xy_range - y0).^2; yv = min(tmp); yi = find((tmp-yv)<1e-6,1,'first');
            yd = (y0 - T.xy_range(yi)) / T.xy_step;
            %{
            [~,xi] = min((T.xy_range - x0).^2);
            xd = (x0 - T.xy_range(xi)) / T.xy_step;
            [~,yi] = min((xy_range - y0).^2);
            yd = (y0 - T.xy_range(yi)) / T.xy_step;
            %}

            rr = yi:T.zoom:size(T.psf_dx,1);
            cc = xi:T.zoom:size(T.psf_dx,2);
            rr = rr(rows(find(rows,1,'first'):end));
            cc = cc(cols(find(cols,1,'first'):end));
        end
        
    end
    
    %% PSF generators
    methods (Access = protected)
       
        function G = genPSF_Analytic(T,xgrid,ygrid,mol)
            global ind;
            
            G = zeros(size(xgrid));
            zi = T.calcIndexDiffZ(mol);
            if zi < 0, return; end;    % `z` out of range
            
            if T.cal.is_biplane
                K = AnalyticPSF.genBiplane(squeeze(xgrid(1,:,:)),squeeze(ygrid(:,1,:)),mol(ind.Z),mol(ind.X),mol(ind.Y),T.cal,T.cal.px);
                for plane=1:size(xgrid,3)
                    [rows,cols] = T.calcIndicesXY(xgrid,ygrid,plane,mol);
                    G(rows,cols,plane) = mol(ind.I) .* squeeze(K.stack(rows,cols,1,plane));
                end
            else
                K = AnalyticPSF.genAstigmatic(xgrid(1,:),ygrid(:,1),mol(ind.Z),mol(ind.X),mol(ind.Y),T.cal,T.cal.px);
                plane = 1;
                [rows,cols] = T.calcIndicesXY(xgrid,ygrid,plane,mol);
                G(rows,cols,plane) = mol(ind.I) .* squeeze(K.stack(rows,cols,1,plane));
            end
        end
        
        function G = genPSF_Spline(T,xgrid,ygrid,mol)
            global ind;
            
            G = zeros(size(xgrid));
            [zi,zd] = T.calcIndexDiffZ(mol);
            if zi < 0, return; end;    % `z` out of range
            zd = (zi + zd) .* ones(size(xgrid));
            
            for plane=1:size(xgrid,3)
                [rows,cols,xd,yd] = T.calcIndicesXY(xgrid,ygrid,plane,mol);
                G(rows,cols,plane) = mol(ind.I) .* interp3(T.psf_bspline{plane},xd,yd,zd(rows,cols,plane),0);
            end
        end
        
        function G = genPSF_Taylor1(T,xgrid,ygrid,mol)
            global ind;
            
            G = zeros(size(xgrid));
            [zi,zd] = T.calcIndexDiffZ(mol);
            if zi < 0, return; end;    % `z` out of range
            
            for plane=1:size(xgrid,3)
                [rows,cols,xd,yd] = T.calcIndicesXY(xgrid,ygrid,plane,mol);
                [xd,yd,rr,cc] = T.calcDiffs(rows,cols,xd(1,1),yd(1,1));
                G(rows,cols,plane) = mol(ind.I) .* (T.psf_dI(rr,cc,zi,plane) + xd.*T.psf_dx(rr,cc,zi,plane) + yd.*T.psf_dy(rr,cc,zi,plane) + zd.*T.psf_dz(rr,cc,zi,plane));
            end
        end
        
        function G = genPSF_Taylor2(T,xgrid,ygrid,mol)
            global ind;
            
            G = zeros(size(xgrid));
            [zi,zd] = T.calcIndexDiffZ(mol);
            if zi < 0, return; end;    % `z` out of range
            
            for plane=1:size(xgrid,3)
                [rows,cols,xd,yd] = T.calcIndicesXY(xgrid,ygrid,plane,mol);
                [xd,yd,rr,cc] = T.calcDiffs(rows,cols,xd(1,1),yd(1,1));
                G(rows,cols,plane) = mol(ind.I) .* (T.psf_dI(rr,cc,zi,plane) + xd.*T.psf_dx(rr,cc,zi,plane) + yd.*T.psf_dy(rr,cc,zi,plane) + zd.*T.psf_dz(rr,cc,zi,plane) ...
                    + (xd*xd/2).*T.psf_dxx(rr,cc,zi,plane) + (yd*yd/2).*T.psf_dyy(rr,cc,zi,plane) + (zd*zd/2).*T.psf_dzz(rr,cc,zi,plane) ...
                    + (xd*yd).*T.psf_dxy(rr,cc,zi,plane) + (xd*zd).*T.psf_dxz(rr,cc,zi,plane) + (yd*zd).*T.psf_dyz(rr,cc,zi,plane));
            end
        end
        
    end
    
    %% Jacobian evaluation
    methods (Access = protected)
        
        function J = calcJacobian_Analytic(T,xgrid,ygrid,mol,tmol)
            imgrid = repmat(reshape(1:size(xgrid,3),[1,1,size(xgrid,3)]),[size(xgrid,1),size(xgrid,2),1]);

            global ind;

            cosfi = cos(T.cal.fi);
            sinfi = sin(T.cal.fi);
            [wx,wy] = AnalyticPSF.defocusGaussian(T.cal,tmol(ind.Z));
            J = zeros(length(xgrid(:)),length(tmol));

            xd = xgrid(:) - tmol(ind.X);
            yd = ygrid(:) - tmol(ind.Y);
            cosfiXd = cosfi .* xd; cosfiYd = cosfi .* yd;
            sinfiYd = sinfi .* yd; sinfiXd = sinfi .* xd;
            first = cosfiXd - sinfiYd; second = sinfiXd + cosfiYd;
            expVal = exp(-0.5.*((first./wx(imgrid(:))).^2 + (second./wy(imgrid(:))).^2));
            oneDivPISS2 = 1 ./ (pi.*wx(imgrid(:)).*wy(imgrid(:)));
            % diff(psf, x0)
            pom1 = first.*cosfi./(wx(imgrid(:)).^2) + second.*sinfi./(wy(imgrid(:)).^2);
            J(:,1) = oneDivPISS2 .* 0.5 .* tmol(ind.I) .* (pom1 .* expVal);
            % diff(psf, y0)
            pom2 = first.*sinfi./(wx(imgrid(:)).^2) + second.*cosfi./(wy(imgrid(:)).^2);
            J(:,2) = oneDivPISS2 .* 0.5 .* tmol(ind.I) .* (pom2 .* expVal);
            % diff(psf, z0)
            pom4 = (tmol(ind.Z) - T.cal.cx(imgrid(:))) ./ (T.cal.d .* (1 + ((tmol(ind.Z) - T.cal.cx(imgrid(:))) ./ T.cal.d).^2)).^2;
            pom5 = (tmol(ind.Z) - T.cal.cy(imgrid(:))) ./ (T.cal.d .* (1 + ((tmol(ind.Z) - T.cal.cy(imgrid(:))) ./ T.cal.d).^2)).^2;
            pom3 = first.^2 .* pom4 + second.^2 .* pom5;
            J(:,3) = oneDivPISS2 .* 0.5 .* tmol(ind.I) .* (expVal .* pom3) ...
                   - oneDivPISS2 .* 0.5 .* tmol(ind.I) .* (expVal .* pom4) ...
                   - oneDivPISS2 .* 0.5 .* tmol(ind.I) .* (expVal .* pom5);
            % diff(psf, I)
            J(:,4) = mol(ind.I) .* expVal .* oneDivPISS2;
            % diff(psf, b)
            J(:,5) = 2 * mol(ind.O);
        end
        
        function J = calcJacobian_Spline(T,xgrid,ygrid,mol,tmol)    % TODO: biplane!
            global ind;
            
            dx = zeros(size(xgrid)); dy = dx; dz = dx; dI = dx;
            [zi,zd] = T.calcIndexDiffZ(mol);
            if zi < 0, return; end;    % `z` out of range
            zd = (zi + zd) .* ones(size(xgrid));
            
            for plane=1:size(xgrid,3)
                [rows,cols,xd,yd] = T.calcIndicesXY(xgrid,ygrid,plane,mol);
                
                dx(rows,cols,plane) = -tmol(ind.I) .* interp3(T.psf_dx_bspline{plane},xd,yd,zd(rows,cols,plane),0);
                dy(rows,cols,plane) = -tmol(ind.I) .* interp3(T.psf_dy_bspline{plane},xd,yd,zd(rows,cols,plane),0);
                dz(rows,cols,plane) =   mol(ind.I) .* interp3(T.psf_dz_bspline{plane},xd,yd,zd(rows,cols,plane),0);
                dI(rows,cols,plane) = 2*mol(ind.I) .* interp3(T.psf_bspline   {plane},xd,yd,zd(rows,cols,plane),0);
                dO = 2*mol(ind.O) .* ones(size(xgrid));
            end
            
            J = [dx(:),dy(:),dz(:),dI(:),dO(:)];
        end
        
        function J = calcJacobian_Taylor1(T,xgrid,ygrid,mol,tmol)
            global ind;
            
            dx = zeros(size(xgrid)); dy = dx; dz = dx; dI = dx;
            [zi,zd] = T.calcIndexDiffZ(mol);
            if zi < 0, return; end;    % `z` out of range
            
            for plane=1:size(xgrid,3)
                [rows,cols,xd,yd] = T.calcIndicesXY(xgrid,ygrid,plane,mol);
                [xd,yd,rr,cc] = T.calcDiffs(rows,cols,xd(1,1),yd(1,1));
                
                dx(rows,cols,plane) = -tmol(ind.I) .* (T.psf_dx(rr,cc,zi,plane) + xd.*T.psf_dxx(rr,cc,zi,plane) + yd.*T.psf_dxy(rr,cc,zi,plane) + zd.*T.psf_dxz(rr,cc,zi,plane));
                dy(rows,cols,plane) = -tmol(ind.I) .* (T.psf_dy(rr,cc,zi,plane) + xd.*T.psf_dyx(rr,cc,zi,plane) + yd.*T.psf_dyy(rr,cc,zi,plane) + zd.*T.psf_dyz(rr,cc,zi,plane));
                dz(rows,cols,plane) =   mol(ind.I) .* (T.psf_dz(rr,cc,zi,plane) + xd.*T.psf_dzx(rr,cc,zi,plane) + yd.*T.psf_dzy(rr,cc,zi,plane) + zd.*T.psf_dzz(rr,cc,zi,plane));
                dI(rows,cols,plane) = 2*mol(ind.I) .* (T.psf_dI(rr,cc,zi,plane) + xd.*T.psf_dx (rr,cc,zi,plane) + yd.*T.psf_dy (rr,cc,zi,plane) + zd.*T.psf_dz (rr,cc,zi,plane));
                dO = 2*mol(ind.O) .* ones(size(xgrid));
            end
            
            J = [dx(:),dy(:),dz(:),dI(:),dO(:)];
        end
        
    end
    
    %% Diffs evaluation
    methods (Access = protected)
        
        function [d1,d2] = calcDataDiffs_P(T,data,estimate)
            d1 = data ./ estimate - 1;
            d2 = data ./ estimate.^2;
        end

        function [d1,d2] = calcDataDiffs_G(T,data,estimate)
            d1 = data - estimate;
            d2 = ones(size(data));
        end
        
    end
    
    %% Public methods
    methods (Access = public)
        
        function T = PSFrefine(psf,low_w,low_h,zoom,z_range,cal,psfEval,jacobianEval,noiseModel)
            % initializing properties
			T.fnGenPsf = str2func(['genPSF_',psfEval]);
			T.fnCalcJacobian = str2func(['calcJacobian_',jacobianEval]);
            T.fnCalcDataDiffs = str2func(['calcDataDiffs_',noiseModel]);
            T.cal = cal;
            T.low_w = low_w;
            T.low_h = low_h;
            T.zoom = zoom;
            T.hi_w = (T.low_w - 1) * T.zoom + 1;
            T.hi_h = (T.low_h - 1) * T.zoom + 1;
            T.hi_z = length(z_range);
            T.z_range = z_range;
            T.z_step = psf.voxel(3);
            if(length(z_range) > 1), T.z_step = abs(z_range(2) - z_range(1)); end;
            T.xy_range = (1:1/T.zoom:T.low_w) - 0.5 - 1/T.zoom/2;
            T.xy_step = 1;
            
            if strcmp(psfEval,'Analytic') & strcmp(jacobianEval,'Analytic'), return; end;
            
            % B-spline approximation of PSF and its partial derivatives
            for zz=1:size(psf,3)
                psf.stack(:,:,zz,:) = psf.stack(:,:,zz,:) ./ sum(sum(sum(psf.stack(:,:,zz,:))));
            end
            spacing = psf.voxel ./ [cal.px,cal.px,T.z_step];
            for ii = 1:size(psf.stack,4)
                T.psf_bspline{ii} = bsarray(psf.stack(:,:,:,ii),'degree',3,'elementSpacing',spacing,'lambda',0);
                %
                T.psf_dx_bspline{ii} = partial(T.psf_bspline{ii},2);
                T.psf_dy_bspline{ii} = partial(T.psf_bspline{ii},1);
                T.psf_dz_bspline{ii} = partial(T.psf_bspline{ii},3);
                
                if strcmp(psfEval,'Taylor1') | strcmp(psfEval,'Taylor2') | strcmp(jacobianEval,'Taylor1')
                    T.psf_dI (:,:,:,ii) = indirectFilter(T.psf_bspline{ii});
                    T.psf_dx (:,:,:,ii) = indirectFilter(T.psf_dx_bspline{ii});
                    T.psf_dy (:,:,:,ii) = indirectFilter(T.psf_dy_bspline{ii});
                    T.psf_dz (:,:,:,ii) = indirectFilter(T.psf_dz_bspline{ii});
                end
                
                if strcmp(psfEval,'Taylor2') | strcmp(jacobianEval,'Taylor1')
                    T.psf_dxx(:,:,:,ii) = indirectFilter(partial(T.psf_dx_bspline{ii},2));
                    T.psf_dxy(:,:,:,ii) = indirectFilter(partial(T.psf_dx_bspline{ii},1));
                    T.psf_dxz(:,:,:,ii) = indirectFilter(partial(T.psf_dx_bspline{ii},3));
                    T.psf_dyx(:,:,:,ii) = indirectFilter(partial(T.psf_dy_bspline{ii},2));
                    T.psf_dyy(:,:,:,ii) = indirectFilter(partial(T.psf_dy_bspline{ii},1));
                    T.psf_dyz(:,:,:,ii) = indirectFilter(partial(T.psf_dy_bspline{ii},3));
                    T.psf_dzx(:,:,:,ii) = indirectFilter(partial(T.psf_dz_bspline{ii},2));
                    T.psf_dzy(:,:,:,ii) = indirectFilter(partial(T.psf_dz_bspline{ii},1));
                    T.psf_dzz(:,:,:,ii) = indirectFilter(partial(T.psf_dz_bspline{ii},3));
                end
            end
        end
        
        function G = genPSF(T,xgrid,ygrid,mol)
            G = T.fnGenPsf(T,xgrid,ygrid,mol);
        end
        
        function J = calcJacobian(T,xgrid,ygrid,mol,tmol)
            J = T.fnCalcJacobian(T,xgrid,ygrid,mol,tmol);
        end
        
        function [d1,d2] = calcDataDiffs(T,data,estimate)
            [d1,d2] = T.fnCalcDataDiffs(T,data,estimate);
        end
        
    end
    
end
