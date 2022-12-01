function [res, resperf, resNoise, fint, fperf, resNoise_perf] = calc_res_val(bval,trifunctie_GT,lb,ub,elements,int_cutoff,nperm,seed,gt_1,gt_2,gt_3,gt_4,gt_5,gt_6)
    D_space = logspace(log10(lb),log10(ub),elements);
    A = exp( -kron( bval',D_space ) );

    fint = zeros(1,nperm);
    fperf = zeros(1,nperm);
    Dstar = zeros(1,nperm);
    D = zeros(1,nperm);
    Dint = zeros(1,nperm);
    resNoise = zeros(1,nperm);
    resNoise_perf = zeros(1,nperm);
    resNoise_D = zeros(1,nperm);
    resNoise_Dstar = zeros(1,nperm);
    resNoise_Dint = zeros(1,nperm);
    
    par = (D_space<1.5e-3);
    int = (D_space<int_cutoff)-par;
    perf = D_space>int_cutoff;          
    rng(seed)
    
    for nNoise = 1:nperm    
        
        gt = [gt_1(nNoise) gt_2(nNoise) gt_3(nNoise) gt_4(nNoise) gt_5(nNoise) gt_6(nNoise)];

        noise = 1.*randn(length(bval),1);

        diff_signal = trifunctie_GT(gt,bval(1:end)') + noise;

        [x] = lsqnonneg(A,[diff_signal]);

        D(nNoise) = sum(x(par).*D_space(par)')/sum(x(par));
        Dint(nNoise) = sum(x(int>0).*D_space(int>0)')/sum(x(int>0));
        Dstar(nNoise) = sum(x(perf).*D_space(perf)')/sum(x(perf));
        
        fint(nNoise) = sum(x'.*int) / sum(x); %if fint(nNoise) == 0; Dint(nNoise) = 0; end;
        fperf(nNoise) = sum(x'.*perf) / sum(x);
        resNoise_perf(nNoise) = abs((fperf(nNoise)*100 - gt(2)*100).^2);
        resNoise_Dstar(nNoise) = abs((Dstar(nNoise)*100 - gt(5)*100).^2);
        resNoise_D(nNoise) = abs((D(nNoise) - gt(3)).^2);
        resNoise_Dint(nNoise) = abs((Dint(nNoise)*100 - gt(4)*100).^2);

        resNoise(nNoise) = abs((fint(nNoise)*100 - gt(1)*100).^2);
    end
    res = mean(resNoise);   
    resperf = mean(resNoise_perf); 
end