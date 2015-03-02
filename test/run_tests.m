function run_tests()
    addpath(genpath('.'));
    addpath(genpath('../src'));
    
    fprintf('\n============ COMMON ============\n');
    fprintf(':: fullAAIinv\n'); fullAAIinv();
    fprintf(':: smallPSF\n'); smallPSF();
    
    fprintf('\n============ BIPLANE ============\n');
    run_all(true);
    
    fprintf('\n\n============ ASTIGMATISM ============\n');
    run_all(false);
end

function run_all(is_biplane)
    fprintf('\n== PSFdeconv ==\n'); run_PSFdeconv(is_biplane);
    fprintf('\n== PSFrefine ==\n'); run_PSFrefine(is_biplane);
end

function run_PSFdeconv(is_biplane)
    fprintf(':: constructor\n'); constructor(is_biplane);
    fprintf(':: convolution\n'); convolution(is_biplane);
    fprintf(':: correlation\n'); correlation(is_biplane);
    fprintf(':: genAAIinv\n'); genAAIinv(is_biplane);
end

function run_PSFrefine(is_biplane)
    fprintf(':: genPSF\n'); genPSF(is_biplane);
    fprintf(':: calcJacobian\n'); calcJacobian(is_biplane);
end