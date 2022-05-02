%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script properties
% name : 
% Description :
% Arguments :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
% 20201711 0.1 Gerald Drenthen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nperm = 10000;
int_cutoff = 4e-3;
lb = .1e-3;
ub = 200e-3;
rng(999)
NSA = 100; % number of voxels, mimicking region-wise fitting
gt_1 = 0 + (0.3 -0) * rand(nperm,1);
gt_2 = 0 + (0.1 -0) * rand(nperm,1);
gt_3 = (1.5e-3-lb)/8 * randn(nperm,1) + (1.5e-3+lb)/2; gt_3(gt_3<lb)=lb; gt_3(gt_3>1.5e-3)=1.5e-3;
gt_4 = (int_cutoff-1.5e-3)/8 * randn(nperm,1) + (int_cutoff+1.5e-3)/2; gt_4(gt_4<1.5e-3)=1.5e-3; gt_4(gt_4>int_cutoff)=int_cutoff;
gt_5 = (ub-int_cutoff)/8 * randn(nperm,1) + (ub+int_cutoff)/2; gt_5(gt_5<int_cutoff)=int_cutoff; gt_5(gt_5>ub)=ub;
S0 = 100;

trifunctie_GT=@(x,bval_tri) ( x(6)* ... 
                         ((1-x(1)-x(2)) * exp(-(bval_tri*x(3))) ...
                        + x(1) * exp(-(bval_tri*(x(3)+x(4)))) ...
                        + x(2) * exp(-(bval_tri*(x(3)+x(4)+x(5))))));
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Define b-value sampling scheme's to assess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bValSet = 0:10:1000;
bval_opt = [0 30 90 210 280 350 580 620 660 680 720 760 980 990 1000];
bval_log = round(logspace(log10(1),log10(101),15));
bval_log(2:6) = bval_log(2:6) + 1;
bval_lin = round(linspace(1,101,15)); 

Aopt = []; Alog = [];
for nn = 1:length(D_space)
    Aopt(:,nn) = exp(-1.*bval_opt.*D_space(nn));
    Alog(:,nn) = exp(-1.*bval_log.*D_space(nn));
    Alin(:,nn) = exp(-1.*bval_lin.*D_space(nn));
end

D_space = logspace(log10(.1e-3),log10(200e-3),200);
par = (D_space<1.5e-3);
int = (D_space<4e-3)-par;
perf = D_space>4e-3;  
rangeComp = par + 2.*int + 3.*perf;

for nNoise = 1:nperm
    rng(nNoise)
    gt = [gt_1(nNoise) gt_2(nNoise) gt_3(nNoise) gt_4(nNoise) gt_5(nNoise) S0];
    finttrue(nNoise) = gt(1);     
    noise = (1/sqrt(NSA)).*randn(length(bValSet),1); 
    
    noise_opt = noise(bval_opt/10+1);
    diff_signal = trifunctie_GT(gt,bval_opt(1:end)') + noise_opt;
    [x] = lsqnonneg(Aopt,diff_signal);
    fintopt(nNoise) = sum(x'.*int) / sum(x);
    
    noise_log = noise(bval_log);
    diff_signal = trifunctie_GT(gt,bval_log(1:end)') + noise_log;
    [x] = lsqnonneg(Alog,diff_signal);
    fintlog(nNoise) = sum(x'.*int) / sum(x);
    
    noise_ref = noise(bval_lin);
    diff_signal = trifunctie_GT(gt,bval_lin(1:end)') + noise_ref;
    [x] = lsqnonneg(Alin,diff_signal);
    fintlin(nNoise) = sum(x'.*int) / sum(x);  
end   

figure; scatter(finttrue,fintopt-finttrue,'.k'); hold on; scatter(finttrue(fintopt==0),fintopt(fintopt==0)-finttrue(fintopt==0),'.','MarkerEdgeColor',[.8 .8 .8]); xlim([0 .3]); ylim([-.3 1])
figure; scatter(finttrue,fintlog-finttrue,'.k'); hold on; scatter(finttrue(fintlog==0),fintlog(fintlog==0)-finttrue(fintlog==0),'.','MarkerEdgeColor',[.8 .8 .8]); xlim([0 .3]); ylim([-.3 1])
figure; scatter(finttrue,fintlin-finttrue,'.k'); hold on; scatter(finttrue(fintlin==0),fintlin(fintlin==0)-finttrue(fintlin==0),'.','MarkerEdgeColor',[.8 .8 .8]); xlim([0 .3]); ylim([-.3 1])
rmse_opt = rms((fintopt-finttrue));
rmse_ref = rms((fintlog-finttrue));
rmse_lin = rms((fintlin-finttrue));

rmse_ref./rmse_opt
rmse_lin./rmse_opt