%analyzing thresholds 
function[mt] = SigmoidThreshold(coeffs, xq, threshold)
    SigmoidFun = GetSigmoid(length(coeffs));
    y_fit = SigmoidFun(coeffs,xq);
    %Converting sigmoid to d'
    dprimeq = norminv(y_fit) - norminv(y_fit(1));
    yq_idx = find(dprimeq >= threshold,1,'first');
    mt = xq(yq_idx);
    if isempty(yq_idx)
        mt = NaN;
    end
end

