function mol = detect_peaks(estI,thresh,xy_zoom,z_range)

    global ind nparams;
    
    z_step = 0;
    if(length(z_range) > 1)
        z_step = abs(z_range(2) - z_range(1));
        z_range = z_range(:);   % always ensure column vector
    end;

    [Ny,Nx,Nz] = size(estI);
    img = zeros(Ny+2,Nx+2,Nz,'single');

    I_mask = ones(3,3,Nz);
    x_mask = repmat([-1 0 1;-1 0 1;-1 0 1],[1,1,Nz]);
    x_mask = x_mask(:);
    y_mask = repmat([-1 -1 -1; 0 0 0; 1 1 1],[1,1,Nz]);
    y_mask = y_mask(:);
    z_mask = zeros(3,3,size(estI,3));
    zind = 1;
    for zval = -floor(Nz/2):+floor(Nz/2)
        z_mask(:,:,zind) = zval .* ones(3,3);
        zind = zind + 1;
    end
    z_mask = z_mask(:);

    img(2:end-1,2:end-1,:) = estI;
    img_flat = sum(img,3);
    img_flat(img_flat<thresh) = 0;
    
    idx = find(img_flat(:)>0);
    [I,J] = ind2sub(size(img_flat),idx);
    nmols = 0;
    intensity = zeros(1000,1);
    xpos = zeros(1000,1); ypos = zeros(1000,1); zpos = zeros(1000,1);
    for ii = 1:length(I)
        xx = J(ii);
        yy = I(ii);
        img_temp = img_flat(yy-1:yy+1,xx-1:xx+1);

        % for local maxima
        if (img_flat(yy,xx) >= max(img_temp(:))) && (img_flat(yy,xx) >= thresh)
            % center of mass
            img_temp = img(yy-1:yy+1,xx-1:xx+1,:);
            img_temp2 = I_mask(:).*img_temp(:);
            %
            nmols = nmols + 1;
            intensity(nmols) = sum(img_temp2);
            xpos(nmols) = sum(x_mask.*img_temp2)./intensity(nmols) + xx - 1;
            ypos(nmols) = sum(y_mask.*img_temp2)./intensity(nmols) + yy - 1;
            zpos(nmols) = sum(z_mask.*img_temp2)./intensity(nmols) + round(Nz/2);
        end
    end
    
    mol = zeros(1,nparams*nmols);
    mol(ind.X:nparams:end) = xpos(1:nmols) ./ xy_zoom;
    mol(ind.Y:nparams:end) = ypos(1:nmols) ./ xy_zoom;
    mol(ind.Z:nparams:end) = z_range(round(zpos(1:nmols))) + z_step*(zpos(1:nmols) - round(zpos(1:nmols)));
    mol(ind.I:nparams:end) = intensity(1:nmols);
    mol(ind.O:nparams:end) = 0.1;

end