classdef GPU
    
    properties(Constant)
        useGPU = true;
    end
    
    methods (Static)
        
        function GA = to(A)
            if GPU.useGPU
                GA = gpuArray(single(A));
            else
                GA = A;
            end
        end
        
        function A = from(GA)
            if GPU.useGPU
                A = gather(GA);
            else
                A = GA;
            end
        end
        
    end
    
end
