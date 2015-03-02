function chi2 = calc_chi2(y,f)
    chi2 = 2.*(f - y - y.*log(f./y));
    chi2 = sum(chi2(:));
end
