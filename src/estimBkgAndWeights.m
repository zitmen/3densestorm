function [bkg,w] = estimBkgAndWeights(A,y,bkgWfactor,bkgSigma,bkgIter)
    if nargin == 2, bkgWfactor = 1.0; end;
    bkg = estimBG(y,bkgSigma,bkgIter);
    eps = min(bkg(:));
    if(eps < 0), bkg = bkg - eps + 0.1; end;
    ww = imresize(bkgWfactor*sqrt(bkg),[A.hi_h,A.hi_w],'bilinear');
    w = zeros(A.hi_h,A.hi_w,A.hi_z);
    for ii = 1:size(y,4)
        w = w + (repmat(ww(:,:,ii),[1,1,A.hi_z,size(y,4)]) ./ size(y,4));
    end
end

function A = estimBG(I,sigma,iterations)
    if(nargin < 2), sigma = 10.0; end;
    if(nargin < 3), iterations = 10; end;
    
    imsiz = [size(I,1),size(I,2)];

    IM = I;
    for i = 1:iterations
        A = zeros(size(IM));
        for n = 1:size(IM,3)
            A(:,:,n) = imgaussian(IM(:,:,n),sigma,imsiz);
        end
        idx = (IM(:) > A(:));
        IM(idx) = A(idx);
    end 
end
