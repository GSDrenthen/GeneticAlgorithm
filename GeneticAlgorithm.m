%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script properties
% name : GeneticAlgorithm.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
% 20201711 0.1 Gerald Drenthen
% 20220305 0.2 Gerald Drenthen
% 20221010 1.0 Gerald Drenthen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clearvars
int_cutoff = 4e-3;
noise_lvl = 1;
S0 = 100;
NbValue = 15;
nperm = 100000;
bValSet = 0:10:1000;
Generations = 250;
PopulationSize = 100;
parents = 20;
offspring = 40;
mutations = 40;
ncores = 12;
allow_multiple_bvals = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Genetic Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gt = [0.30 0.03 8e-4 2e-3 5e-2 S0];

trifunctie_GT=@(x,bval_tri) ( x(6)* ... 
                         ((1-x(1)-x(2)) * exp(-(bval_tri*x(3))) ...
                        + x(1) * exp(-(bval_tri*(x(3)+x(4)))) ...
                        + x(2) * exp(-(bval_tri*(x(3)+x(4)+x(5))))));

for random = 1:50
    rng(random);

    Population = zeros(length(bValSet),PopulationSize);
    warning('off')

    Evolution = zeros(Generations,1);
    Fitness = zeros(PopulationSize,1);

    no_improvement = 0;

    % Generate population
    for nn = 1:PopulationSize
        indx = randi(length(bValSet),[1 NbValue]); while nnz(unique(indx)) ~= NbValue; indx = randi(length(bValSet),[1 NbValue]); end
        if nnz(indx==1) == 0; indx(1) = 1; end
        Population(indx,nn) = 1;
    end

    lb = .1e-3;
    ub = 200e-3;
    gt_1 = 0 + (0.3 -0) * rand(nperm,1);
    gt_2 = 0 + (0.1 -0) * rand(nperm,1);
    gt_3 = (1.5e-3-lb)/8 * randn(nperm,1) + (1.5e-3+lb)/2; gt_3(gt_3<lb)=lb; gt_3(gt_3>1.5e-3)=1.5e-3;
    gt_4 = (int_cutoff-1.5e-3)/8 * randn(nperm,1) + (int_cutoff+1.5e-3)/2; gt_4(gt_4<1.5e-3)=1.5e-3; gt_4(gt_4>int_cutoff)=int_cutoff;
    gt_5 = (ub-int_cutoff)/8 * randn(nperm,1) + (ub+int_cutoff)/2; gt_5(gt_5<int_cutoff)=int_cutoff; gt_5(gt_5>ub)=ub;

    delete(gcp('nocreate'))
    parpool(ncores);        

    for nEvol = 1:Generations
        % Determine fitness
        if nnz(Population(1,:))<100
            for popsize = 1:100
                if Population(1,popsize) == 0
                    indx = find(Population(:,popsize));
                    Population(indx(randi(length(indx),[1 1])),popsize) = 0;
                    Population(1,popsize) = 1;
                end
            end
        end

        pctRunOnAll warning('off')

        parfor nn = 1:PopulationSize
            [Fitness(nn), resnoise(nn,:), fint(nn,:)] = calc_res(bValSet(Population(:,nn)>0),trifunctie_GT,S0,noise_lvl,lb,ub,200,int_cutoff,nperm,nEvol,gt_1,gt_2,gt_3,gt_4,gt_5); 
        end

        [~, indx] = sort(Fitness,'asc');
        Evolution(nEvol) = Fitness(indx(1));
        best_method = indx(1);

        if nEvol == Generations
            error_fint = resnoise(indx(1),:);
            meas_fint = fint(indx(1),:);
        end

        if indx(1) == 1
            no_improvement = no_improvement + 1;
        else
            no_improvement = 0;
        end
        
        disp(['Evolution: ' num2str(nEvol) ', Best Fitness is: ' num2str(Fitness(indx(1))) ', Targeted mutations'])
        
        % Select good parents 
        Population(:,1:parents) = Population(:,indx(1:parents));
        Population(:,parents+1:PopulationSize) = zeros(length(bValSet),PopulationSize-parents);

        % New genes
        if no_improvement == 10
            disp('Introducing new genetics')
            Population(:,(parents+1):PopulationSize) = zeros(length(bValSet),PopulationSize-parents);
            for nn = (parents+1):PopulationSize
                inx = randi(length(bValSet),[1 NbValue]); while nnz(unique(inx)) ~= NbValue; inx = randi(length(bValSet),[1 NbValue]); end
                if nnz(indx==1) == 0; indx(1) = 1; end
                Population(:,nn) = zeros(length(bValSet),1);
                Population(inx,nn) = 1;
            end    
            no_improvement = 0;
        else    

            % Create offspring 
            for nn = 1:offspring
                nPop = randi(parents,[2 1]); while nPop(1) == nPop(2); nPop = randi(parents,[2 1]); end

                Child = zeros(length(bValSet),1);
                if allow_multiple_bvals == 1
                    indx1 = []; for i = 0:max(Population(:,nPop(1)))-1; indx1 = [indx1; find(Population(:,nPop(1))>i)]; end                
                    indx2 = []; for i = 0:max(Population(:,nPop(2)))-1; indx2 = [indx2; find(Population(:,nPop(2))>i)]; end                
                    indx1 = sort(indx1);
                    indx2 = sort(indx2);
                else                
                    indx1 = find(Population(:,nPop(1)));
                    indx2 = find(Population(:,nPop(2)));
                end

                for ll = 1:NbValue
                      if randi([0 1],1) == 0
                           Child(indx1(ll)) = 1;
                      else
                           Child(indx2(ll)) = 1;
                      end
                end

                %if randi([0 1],1) == 0
                    numMut = randi(3,1);
                    for ll = 1:numMut
                         indx = find(Child(2:end));
                         Child(1+indx(randi(length(indx),1))) = 0;
                         Child(randi(length(bValSet),1)) = 1;
                         if allow_multiple_bvals == 1
                            while sum(Child) > NbValue; tmp=1+randi(length(bValSet),1); Child(tmp) = Child(tmp) -1; end
                            while sum(Child) < NbValue; tmp=randi(length(bValSet),1); Child(tmp) = Child(tmp) + 1; end
                         else
                            while nnz(Child) > NbValue; Child(1+indx(randi(length(indx),1))) = 0; end
                            while nnz(Child) < NbValue; Child((randi(length(bValSet),1))) = 1; end
                         end
                    end
                %end
                Population(Child>0,parents+nn) = 1;

            end


            for nn = 1:mutations
                if nn < floor(mutations/2) % Random Mutations 
                    nPop = randi(parents+offspring,[1 1]);
                    Mutation = Population(:,nPop(1))>0;
                    numMut = randi(NbValue-1,1);
                    if allow_multiple_bvals == 1
                        indx = []; for i = 0:max(Mutation)-1; indx = [indx; find(Mutation>i)]; end                
                        indx = sort(indx(2:end));
                        tmp=indx(randi(length(indx),numMut,1)); Mutation(tmp) = Mutation(tmp) - 1;
                        while sum(Mutation) > NbValue; tmp=find(Mutation); indx=tmp(randi(nnz(Mutation))); Mutation(indx) = Mutation(indx) - 1; end
                        while sum(Mutation) < NbValue; tmp=randi(length(bValSet),1); Mutation(tmp) = Mutation(tmp) + 1; end                            
                    else
                        indx = find(Mutation(2:end));                        
                        Mutation(indx+1) = 0;
                        Mutation(randi(length(bValSet),numMut,1)) = 1;
                        while nnz(Mutation) > NbValue; Mutation(randi(length(bValSet),1)) = 0; end
                        while nnz(Mutation) < NbValue; Mutation(randi(length(bValSet),1)) = 1; end
                    end
                else % Targeted mutations
                    nPop = randi(parents+offspring,[1 1]);
                    Mutation = Population(:,nPop(1))>0;
                    numMut = randi(NbValue-1,1);
                    if allow_multiple_bvals == 1
                        indx = []; for i = 0:max(Mutation)-1; indx = [indx; find(Mutation>i)]; end                
                        indx = sort(indx(2:end));                        
                    else
                        indx = find(Mutation(2:end));
                    end                    
                    TargMut = randi(4,numMut,1);
                    i = 1;
                    for ll = 1:length(TargMut)
                        if randi(2,1) == 1
                            TargMut(ll) = TargMut(ll)*-1;
                        end
                    end
                    indx_mut = zeros(NbValue-1,1);
                    for ll = 1:length(TargMut)
                        indx_mut(randi(NbValue-1,1)) = TargMut(ll);
                    end
                    indx_mut = indx_mut + indx;
                    for ll = 1:NbValue-1   
                        if indx_mut(ll) > length(bValSet)-1; indx_mut(ll) = length(bValSet)-1; end
                        if indx_mut(ll) < 1; indx_mut(ll) = 1; end
                    end
                    if allow_multiple_bvals == 1
                       if max(Mutation(indx))>1
                          for i = 1:length(indx);  Mutation(indx(i)) = Mutation(indx(i)) - 1; end
                       else
                           Mutation(indx) = Mutation(indx) - 1;
                       end
                       Mutation(indx_mut) = Mutation(indx_mut) + 1;
                       while sum(Mutation) > NbValue; tmp=find(Mutation); indx=tmp(randi(nnz(Mutation))); Mutation(indx) = Mutation(indx) - 1; end
                       while sum(Mutation) < NbValue; tmp=randi(length(bValSet),1); Mutation(tmp) = Mutation(tmp) + 1; end
                    else
                        Mutation(indx+1) = 0;
                        Mutation(indx_mut+1) = 1;                        
                        while nnz(Mutation) > NbValue; Mutation(randi(length(bValSet),1)) = 0; end
                        while nnz(Mutation) < NbValue; Mutation(randi(length(bValSet),1)) = 1; end
                    end            
                end
                Population(Mutation,parents+offspring+nn) = 1;
            end    
        end
        Population_3D(:,:,nEvol) = Population;
    end
    save('GeneticAlgorithm.mat')
end