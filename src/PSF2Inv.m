classdef PSF2Inv
    
    properties (SetAccess = private)
        A;
    end
    
    properties (SetAccess = public)
        fpsfi;
        hi_h; hi_w;
        fnHessMult;
    end
    
    %% Private helpers
    methods (Static, Access = private)
    
        function M = block_mul(A,B,fn)
        % A...represents a PSF
        % B...represents a PSF when building AAIinv and an image otherwise
            [ma,na] = size(A);
            [mb,nb] = size(B);
            if mb ~= na, error('dimensions dont agree'); end;
            M = cell(ma,nb);
            for i = 1:ma
                for j = 1:nb
                    M{i,j} = zeros(size(B{i,j}));
                    for k = 1:na
                        M{i,j} = M{i,j} + fn(A{i,k}, B{k,j});
                    end
                end
            end
        end
        
    end
    
    %% Block matrix operations
    methods (Static, Access = public)

        function [U,D,L] = SchurDecomposition(AtA)
            n = size(AtA,1);
            M = cell(n,n);  % compact form of U,D,L; by defitnition it can fit all into a single matrix
            U = cell(n,n); D = cell(n,n); L = cell(n,n);
            for r = n:-1:1
                for c = n:-1:1
                    % calculate the general recursive formula
                    k = max(r,c);
                    M{r,c} = AtA{r,c};
                    for s=n:-1:k+1
                        M{r,c} = M{r,c} - M{r,s}.*M{s,c}./M{s,s};
                    end
                    % fill U,L,D respectively
                    if r == c
                        D{r,c} = M{r,c};
                        L{r,c} = ones(size(AtA{r,c}));  % fft of identity matrix
                        U{r,c} = ones(size(AtA{r,c}));  % fft of identity matrix
                    elseif r > c
                        D{r,c} = zeros(size(AtA{r,c}));
                        L{r,c} = M{r,c}./M{k,k};
                        U{r,c} = zeros(size(AtA{r,c}));
                    else % r < c
                        D{r,c} = zeros(size(AtA{r,c}));
                        L{r,c} = zeros(size(AtA{r,c}));
                        U{r,c} = M{r,c}./M{k,k};
                    end
                end
            end
        end
        
        function M = BlockTranspose(A)
            [m,n] = size(A);
            M = cell(m,n);
            for i = 1:m
                for j = 1:n
                    M{j,i} = A{i,j}';
                end
            end
        end
        
        function M = MatrixBlockMultilply(A,B)
            M = PSF2Inv.block_mul(A,B,@(x,y) x.*y);
        end

        function Dinv = DiagonalBlockInverse(D)
            m = length(D);
            Dinv = cell(m,m);
            for i=1:m
                for j=1:m
                    if i==j
                        Dinv{i,j} = 1./D{i,j};
                    else
                        Dinv{i,j} = zeros(size(D{i,j}));
                    end
                end
            end
        end

        function Linv = LowerTriangularBlockInverse(L)
            [m,n] = size(L);
            % prepare the inverse matrix
            Linv = cell(m,n);
            for i = 1:m
                for j = 1:n
                    if i == j
                        Linv{i,j} = ones(size(L{i,j}));  % fft of identity matrix
                    else
                        Linv{i,j} = zeros(size(L{i,j}));
                    end
                end
            end
            % forward
            for j=1:m-1
                for i=j+1:m
                    tmp = L{i,j};
                    for k=1:n
                        L{i,k} = L{i,k} - tmp.*L{j,k};
                        Linv{i,k} = Linv{i,k} - tmp.*Linv{j,k};
                    end
                end
            end
        end

        function Uinv = UpperTriangularBlockInverse(U)
            [m,n] = size(U);
            % prepare the inverse matrix
            Uinv = cell(m,n);
            for i = 1:m
                for j = 1:n
                    if i == j
                        Uinv{i,j} = ones(size(U{i,j}));  % fft of identity matrix
                    else
                        Uinv{i,j} = zeros(size(U{i,j}));
                    end
                end
            end
            % backward
            for j=m:-1:2
                for i=j-1:-1:1
                    tmp = U{i,j};
                    for k=1:n
                        U{i,k} = U{i,k} - tmp.*U{j,k};
                        Uinv{i,k} = Uinv{i,k} - tmp.*Uinv{j,k};
                    end
                end
            end
        end
        
    end
    
    %% Single block multiplication methods
    methods (Access = protected)
        
        function M = mul_fft(T,U)
            UU = cell(size(U,3),1);
            for z=1:size(U,3)
                UU{z} = fft2(U(:,:,z));
            end
            MM = T.MatrixBlockMultilply(T.fpsfi,UU);
            M = GPU.to(zeros(size(U)));
            for z=1:size(U,3)
                M(:,:,z) = real(ifft2(MM{z}));
            end
        end
        
        function M = mul_core(T,U)
            UU = cell(size(U,3),1);
            for z=1:size(U,3)
                UU{z} = U(:,:,z);
            end
            MM = T.block_mul(T.fpsfi,UU,@(h,x) conv2(x,h,'same'));
            M = GPU.to(zeros(size(U)));
            for z=1:size(U,3)
                M(:,:,z) = MM{z};
            end
        end
        
        function M = mul_zero(T,U)
            M = T.mul_core(U);
        end
        
        function M = mul_symmetric(T,U)
            M = T.A.unpaddZoomedImg(T.mul_core(T.A.paddZoomedImg(U,'symmetric')));
        end
        
        function M = mul_circular(T,U)
            M = T.A.unpaddZoomedImg(T.mul_core(T.A.paddZoomedImg(U,'circular')));
        end
        
        function M = mul_replicate(T,U)
            M = T.A.unpaddZoomedImg(T.mul_core(T.A.paddZoomedImg(U,'replicate')));
        end
        
    end
    
    %% Public methods
    methods (Access = public)
        
        function T = PSF2Inv(A,AtA,hessMult)
            T.A = A;
            % init
            T.hi_h = A.hi_h;
            T.hi_w = A.hi_w;
            % calculate the inverse
            [U,D,L] = T.SchurDecomposition(T.BlockTranspose(AtA));
            Ui = T.UpperTriangularBlockInverse(U);
            Di = T.DiagonalBlockInverse(D);
            Li = T.LowerTriangularBlockInverse(L);
            T.fpsfi = T.MatrixBlockMultilply(T.MatrixBlockMultilply(Li,Di),Ui);
            % convert back from the Fourier domain, if needed for `hessMult`
            T.fnHessMult = str2func(['mul_',hessMult]);
            if ~strcmp(hessMult,'fft')
                for i=1:size(T.fpsfi,1)
                    for j=1:size(T.fpsfi,2)
                        T.fpsfi{i,j} = real(ifftshift(ifft2(T.fpsfi{i,j})));
                    end
                end
            end
        end
        
        function M = mtimes(T,U)
            if(isa(U,'PSF2Inv'))
                error('`PSF2Inv * PSF2Inv` is not supported at the moment!');
            elseif((size(U,1) == T.hi_h) & (size(U,2) == T.hi_w))
                M = T.fnHessMult(T,U);
            else
                error('Length of the vector doesn''t match either # of columns nor # of rows of the matrix!');
            end
        end
        
    end
end
