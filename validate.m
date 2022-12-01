clearvars
close all

int_cutoff = 4e-3;
trifunctie_GT=@(x,bval_tri) ( x(6)* ... 
                         ((1-x(1)-x(2)) * exp(-(bval_tri*x(3))) ...
                        + x(1) * exp(-(bval_tri*(x(3)+x(4)))) ...
                        + x(2) * exp(-(bval_tri*(x(3)+x(4)+x(5))))));
bval_opt = [0 30 90 210 280 350 580 620 660 680 720 760 980 990 1000];

bValSet = 0:10:1000;

bval_log = round(logspace(log10(1),log10(101),15));
bval_log(2:6) = bval_log(2:6) + 1;
bval_lin = round(linspace(1,101,15)); 
bval_lit = [0 20 60 150 160 170 190 200 260 440 560 600 700 980 1000];
lb = .1e-3;
ub = 200e-3;

niter_SNR = 100; % SNR is 10 times higher (1/sqrt(N))
niter = 1000;

rng(1)

gt_1_list = 0 + (0.3-0) .*rand(niter,1);
i = 1;
for SNR = [25 50 100]
    
    for N = 1:niter
        rng(N)
        gt_1 = gt_1_list(N).*ones(niter_SNR,1);
        gt_2 = 0 + (0.1 -0) * rand(niter_SNR,1);
        gt_3 = (1.5e-3-lb)/8 * randn(niter_SNR,1) + (1.5e-3+lb)/2; gt_3(gt_3<lb)=lb; gt_3(gt_3>1.5e-3)=1.5e-3;
        gt_4 = (int_cutoff-1.5e-3)/8 * randn(niter_SNR,1) + (int_cutoff+1.5e-3)/2; gt_4(gt_4<1.5e-3)=1.5e-3; gt_4(gt_4>int_cutoff)=int_cutoff;
        gt_5 = (ub-int_cutoff)/8 * randn(niter_SNR,1) + (ub+int_cutoff)/2; gt_5(gt_5<int_cutoff)=int_cutoff; gt_5(gt_5>ub)=ub;
        gt_6 = SNR.*ones(niter_SNR,1);
    
        [res, resperf, resNoise, fint, fperf, resNoise_perf] = calc_res_val(bval_opt,trifunctie_GT,lb,ub,200,int_cutoff,niter_SNR,N,gt_1,gt_2,gt_3,gt_4,gt_5,gt_6);
        fint_opt(N) = mean(fint);
        bias_opt(N) = (mean(fint)-mean(gt_1));
        bias_opt_perf(N) = (mean(fperf)-mean(gt_2));
        
        [res, resperf, resNoise, fint, fperf, resNoise_perf] = calc_res_val(bValSet(bval_log),trifunctie_GT,lb,ub,200,int_cutoff,niter_SNR,N,gt_1,gt_2,gt_3,gt_4,gt_5,gt_6);
        fint_log(N) = mean(fint);
    
        bias_log(N) = (mean(fint)-mean(gt_1));
        bias_log_perf(N) = (mean(fperf)-mean(gt_2));
    
        [res, resperf, resNoise, fint, fperf, resNoise_perf] = calc_res_val(bValSet(bval_lin),trifunctie_GT,lb,ub,200,int_cutoff,niter_SNR,N,gt_1,gt_2,gt_3,gt_4,gt_5,gt_6);
        fint_lin(N) = mean(fint);
    
        bias_lin(N) = (mean(fint)-mean(gt_1));    
        bias_lin_perf(N) = (mean(fperf)-mean(gt_2));
    
        [res, resperf, resNoise, fint, fperf, resNoise_perf] = calc_res_val(bval_lit,trifunctie_GT,lb,ub,200,int_cutoff,niter_SNR,N,gt_1,gt_2,gt_3,gt_4,gt_5,gt_6);
        fint_lit(N) = mean(fint);
        
        bias_lit(N) = (mean(fint)-mean(gt_1));     
        bias_lit_perf(N) = (mean(fperf)-mean(gt_2));
    end

    RMSElog(i) = sqrt(mean(bias_log.^2));
    RMSElin(i) = sqrt(mean(bias_lin.^2));
    RMSElit(i) = sqrt(mean(bias_lit.^2));
    RMSEopt(i) = sqrt(mean(bias_opt.^2));
    
    RMSElog_p(i) = sqrt(mean(bias_log_perf.^2));
    RMSElin_p(i) = sqrt(mean(bias_lin_perf.^2));
    RMSElit_p(i) = sqrt(mean(bias_lit_perf.^2));
    RMSEopt_p(i) = sqrt(mean(bias_opt_perf.^2));
   
    i = i + 1;

end

RMSElog./RMSEopt
RMSElin./RMSEopt
RMSElit./RMSEopt