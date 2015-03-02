function mollist = refine(cfg,A,y,mol)

    global ind nparams;
    
    % get initial localizations
    if isempty(mol), mollist = []; return; end;
    
    % run the Maximum-Likelihood Fitting using Levenberg-Marquardt algorithm
    mol = transformInverse(mol);
    [mol,gen,~,iter] = lm_fitting(cfg,A,mol,y);
    mol = transform(mol);
    
    % throw away molecules with out-of-range parameters
    xpos = mol(ind.X:nparams:end);
    ypos = mol(ind.Y:nparams:end);
    zpos = mol(ind.Z:nparams:end);
    molI = mol(ind.I:nparams:end);
    molO = mol(ind.O:nparams:end);
    zz = zeros(1,length(zpos));
    for zi=1:length(zpos)
        [~,zz(zi)] = min(abs(zpos(zi) - A.z_range));
    end
    yy = round(ypos); xx = round(xpos);
    in_y = (yy > 0) & (yy <= size(y,1));
    in_x = (xx > 0) & (xx <= size(y,2));
    z_step = A.z_step; if z_step == 0, z_step = 100; end; % allow only a limited range of sigmas in 2D case
    dzz = (zpos - A.z_range(zz)) ./ z_step;
    in_z = (dzz >= -1) & (dzz <= +1);
    in_xyz = in_x & in_y & in_z;
    
    if (sum(in_xyz) == 0), mollist = []; return; end;
    
    % fill the results
    mollist.X = xpos(in_xyz) .* A.cal.px;   % px --> nm
    mollist.Y = ypos(in_xyz) .* A.cal.px;   % px --> nm
    mollist.Z = zpos(in_xyz);   % nm
    mollist.I = molI(in_xyz);   % photons
    mollist.O = molO(in_xyz);   % photons
    %mollist.O = molO(in_xyz); + offsetFromBkg(mollist.X,mollist.Y,A.cal.px,estimBkgAndWeights(A,GPU.from(y-gen),0,cfg.bkg.blur_sigma,cfg.bkg.update_iters));   % photons
    
    fprintf('Refinement: %d iterations (%d removed, %d kept)\n',iter,sum(~in_xyz),sum(in_xyz));
    
end

function O = offsetFromBkg(X,Y,pixelsize,bkg)
    xx = round(X ./ pixelsize); xx(xx < 1) = 1; xx(xx > size(bkg,2)) = size(bkg,2);
    yy = round(Y ./ pixelsize); yy(yy < 1) = 1; yy(yy > size(bkg,1)) = size(bkg,1);
    idx = sub2ind(size(bkg),yy,xx);
    O = bkg(idx);
end

function tmol = transform(mol)
    global ind nparams;
    tmol = mol;
    tmol(ind.I:nparams:end) = mol(ind.I:nparams:end).^2;
    tmol(ind.O:nparams:end) = mol(ind.O:nparams:end).^2;
end

function mol = transformInverse(tmol)
    global ind nparams;
    mol = tmol;
    mol(ind.I:nparams:end) = sqrt(abs(tmol(ind.I:nparams:end)));
    mol(ind.O:nparams:end) = sqrt(abs(tmol(ind.O:nparams:end)));
end

function [mol,gen,chi2,iter] = lm_fitting(cfg,A,mol,img)

    global ind nparams;
    [xgrid,ygrid] = meshgrid((1:size(img,2))-0.5,(1:size(img,1))-0.5);
    if cfg.calibration.is_biplane
        xgrid = repmat(xgrid,[1,1,2]); ygrid = repmat(ygrid,[1,1,2]);
    end
    nmols = length(mol) / nparams;
    img = GPU.from(img);
    gen = calcG(A,xgrid,ygrid,mol);
    chi2 = 1e12;
    chi2_prev = Inf;
    bkg = estimBkgAndWeights(A,img-gen,0,cfg.bkg.blur_sigma,cfg.bkg.update_iters);
    bkgc = initBkgCounts(xgrid,ygrid,mol,cfg.fitradius);
    
    iter = 0;
    total = 1;
    lambdas = cfg.lambda.*ones(1,nmols);
    active = ones(1,nmols);
    while (iter < cfg.max_iter) & (total > 0) & (abs(chi2-chi2_prev) > cfg.min_improvement)
        iter = iter + 1;
        
        if (mod(iter,cfg.bkg_update) == 0)
            bkg = estimBkgAndWeights(A,img-gen,0,cfg.bkg.blur_sigma,cfg.bkg.update_iters);
        end
        
        chi2 = calc_chi2(img,gen+bkg);
        chi2_prev = chi2;
        
        total = 0;
        active_ind = find(active);
        for mi = active_ind
            cmol = mol((mi-1)*nparams+1:mi*nparams);
            [umol,updated,gen,bkgc,uchi2] = update(A,xgrid,ygrid,img,gen,bkg,bkgc,cfg.fitradius,cmol,chi2,lambdas(mi));
            if updated
                mol((mi-1)*nparams+1:mi*nparams) = umol;
                if (abs(uchi2 - chi2) / uchi2) < cfg.min_mol_update   % converged?
                    % also if our of bounds, chi2==uchi2 --> converged
                    active(mi) = 0;
                end
                if(uchi2 < chi2)
                    if lambdas(mi) > cfg.lambda
                        lambdas(mi) = lambdas(mi) / cfg.lambda_factor;
                    end
                else
                    lambdas(mi) = lambdas(mi) * cfg.lambda_factor;
                    if lambdas(mi) > cfg.max_lambda
                        active(mi) = 0;
                    end
                end
                chi2 = uchi2;
                total = total + 1;
            end
        end
        
        if((chi2_prev - chi2) < cfg.max_improvement_full_step)
            active(:) = 0;
        end
        
        if((cfg.debug > 0) && (mod(iter,cfg.debug)==0))
            tmol = transform(mol);
            if size(gen,3) > 1
                tmp_gen = [gen(:,:,1),gen(:,:,2)];
                tmp_bkg = [bkg(:,:,1),bkg(:,:,2)];
            else
                tmp_gen = gen;
                tmp_bkg = bkg;
            end
            figure(1),imagesc(tmp_gen+tmp_bkg),title('total estimate'),colorbar,drawnow;
            figure(2),imagesc(tmp_bkg),title('background estimate'),colorbar,drawnow;
            fprintf('==== iter %d ====\n',iter);
            fprintf('CHI2=%f\n',chi2);
            for mi = 1:nmols
                cmol = tmol((mi-1)*nparams+1:mi*nparams);
                fprintf('Molecule %d\n',mi);
                fprintf('x=%f\n',cmol(ind.X));
                fprintf('y=%f\n',cmol(ind.Y));
                fprintf('z=%f\n',cmol(ind.Z));
                fprintf('I=%f\n',cmol(ind.I));
                fprintf('off=%f\n',cmol(ind.O));
                fprintf('\n');
            end
        end
        
    end
    
    if(cfg.debug > 0)
        tmol = transform(mol);
        if size(gen,3) > 1
            tmp_gen = [gen(:,:,1),gen(:,:,2)];
            tmp_bkg = [bkg(:,:,1),bkg(:,:,2)];
        else
            tmp_gen = gen;
            tmp_bkg = bkg;
        end
        figure(1),imagesc(tmp_gen+tmp_bkg),title('total estimate'),colorbar,drawnow;
        figure(2),imagesc(tmp_bkg),title('background estimate'),colorbar,drawnow;
        fprintf('==== FINAL ====\n');
        fprintf('CHI2=%f\n',chi2);
        for mi = 1:nmols
            cmol = tmol((mi-1)*nparams+1:mi*nparams);
            fprintf('Molecule %d\n',mi);
            fprintf('x=%f\n',cmol(ind.X));
            fprintf('y=%f\n',cmol(ind.Y));
            fprintf('z=%f\n',cmol(ind.Z));
            fprintf('I=%f\n',cmol(ind.I));
            fprintf('off=%f\n',cmol(ind.O));
            fprintf('\n');
        end
    end
    fprintf('Final: l2err = %f\n',norm(img(:) - gen(:) - bkg(:)));

end

function [umol,updated,ugen,ubkgc,uchi2] = update(A,xgrid,ygrid,img,gen,bkg,bkgc,fitradius,mol,chi2,lambda)
% prefix u = updated
%
% img - raw image
% gen - image generated from the estimated parameters
% bkg - estimated background
% bkgc - count of fitting windows per position
% mol - parameters for a single molecule
% chi2 - chi2 given the current parameters in `mol`
    
    global clamp ind nparams;
    clamps = transformInverse(clamp);
    
    % prepare the fitting region
    [rows,cols] = getRegion(xgrid,ygrid,mol(ind.X),mol(ind.Y),fitradius);
    grid_x = xgrid(rows,cols,:); grid_y = ygrid(rows,cols,:);
    
    img_reg = img(rows,cols,:);
    gen_bkg = gen(rows,cols,:) + bkg(rows,cols,:);
    
    % calculate the Jacobian and Hessian to get the update vector
    J = A.calcJacobian(grid_x,grid_y,mol,transform(mol));
    
    % correct the derivative by the bkgc factor
    dbkgc = ones(size(grid_x)) ./ max(bkgc(rows,cols,:),1);
    J(:,ind.O) = J(:,ind.O) .* dbkgc(:);
    
    [Jtmp,Htmp] = A.calcDataDiffs(img_reg, gen_bkg);
    H = (J.*repmat(Htmp(:),[1,length(mol)]))'*J;
    
    L = ones(size(H));
    L(1:nparams+1:end) = (1 + lambda);  % L-M diagonal
    
    alpha = H .* L;
    beta = J'*Jtmp(:);
    delta = (alpha \ beta)';
    
    % update
    umol = mol + delta ./ (1 + abs(delta) ./ clamps);
    
    % check if out-of-bounds
    bx = (umol(ind.X) > 0 & umol(ind.X) < A.hi_w);
    by = (umol(ind.Y) > 0 & umol(ind.Y) < A.hi_h);
    bz = (umol(ind.Z) > (min(A.z_range) - A.z_step) & umol(ind.Z) < (max(A.z_range) + A.z_step));
    
    ugen = gen;
    ubkgc = bkgc;
    if bx & by & bz
        % mol(ind.O)^2 --> transform()
        ugen(rows,cols,:) = ugen(rows,cols,:) - calcG(A,grid_x,grid_y,mol) - mol(ind.O)^2 ./ max(bkgc(rows,cols,:),1);
        ubkgc(rows,cols,:) = ubkgc(rows,cols,:) - 1;
        
        % TODO: this could be more efficient if we only check for
        %       [x,y] going  over/under .5
        [rows,cols] = getRegion(xgrid,ygrid,umol(ind.X),umol(ind.Y),fitradius);
        
        ubkgc(rows,cols,:) = ubkgc(rows,cols,:) + 1;
        ugen(rows,cols,:) = ugen(rows,cols,:) + calcG(A,xgrid(rows,cols,:),ygrid(rows,cols,:),umol) + umol(ind.O)^2 ./ max(ubkgc(rows,cols,:),1);
    end
    uchi2 = calc_chi2(img,ugen+bkg);   % if out of bounds, chi2==uchi2

    updated = 1;
    if(isnan(uchi2) || (uchi2 > chi2))
        umol = mol;
        uchi2 = chi2;
        ugen = gen;
        ubkgc = bkgc;
        updated = 0;    % not updated
    end

end

function [rows,cols] = getRegion(xgrid,ygrid,xpos,ypos,fitradius)
    [~,xi] = min((xgrid(1,:,1) - xpos).^2);
    [~,yi] = min((ygrid(:,1,1) - ypos).^2);
    cols = max(xi-fitradius,1):min(xi+fitradius,size(xgrid,2));
    rows = max(yi-fitradius,1):min(yi+fitradius,size(xgrid,1));
end

function bkgc = initBkgCounts(xgrid,ygrid,mol,fitradius)
    global ind nparams;
    bkgc = zeros(size(xgrid));
    for mi = 1:(length(mol)/nparams)
        [rows,cols] = getRegion(xgrid,ygrid,mol((mi-1)*nparams+ind.X),mol((mi-1)*nparams+ind.Y),fitradius);
        bkgc(rows,cols,:) = bkgc(rows,cols,:) + 1;
    end
end

function G = calcG(A,xgrid,ygrid,mol)
    global nparams;
    tmol = transform(mol);
    G = zeros(size(xgrid));
    nmols = length(tmol)/nparams;
    for mi = 1:nmols
        G = G + A.genPSF(xgrid,ygrid,tmol((mi-1)*nparams+1:mi*nparams));
    end
end
