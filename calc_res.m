function [res, resNoise, fint] = calc_res(bval,trifunctie_GT,S0,noise_lvl,lb,ub,elements,int_cutoff,nperm,seed,gt_1,gt_2,gt_3,gt_4,gt_5)
    D_space = logspace(log10(lb),log10(ub),elements);
    A = [];
    for nn = 1:length(D_space)
        A(:,nn) = exp(-1.*[bval].*D_space(nn));
    end

    par = (D_space<1.5e-3);
    int = (D_space<int_cutoff)-par;
    perf = D_space>int_cutoff;      
    rangeComp = par + 2.*int + 3.*perf;

    rng(seed)
    
    for nNoise = 1:nperm       
        gt = [gt_1(nNoise) gt_2(nNoise) gt_3(nNoise) gt_4(nNoise) gt_5(nNoise) S0];
        noise = noise_lvl.*randn(length(bval),1);
        diff_signal = trifunctie_GT(gt,bval(1:end)') + noise;
        [x,resnorm,resid,exitflag,output,lambda] = lsqnonneg(A,diff_signal);
        fint(nNoise) = sum(x'.*int) / sum(x);
        resNoise(nNoise) = abs((fint(nNoise)*100 - gt(1)*100).^2);
    end
    res = mean(resNoise);   
end