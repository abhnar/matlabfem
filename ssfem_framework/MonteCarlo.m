classdef MonteCarlo< SFEM 
    properties
        
        
       
       
        epsilon_r;
        
       
        
       
        
    end
    

    methods
        
        
        function assignEffectiveMaterialVariation(mc, kle, level)
            if(nargin ==2) level = 'full';  end;

            assignEffectiveMaterialVariation@SFEM(mc, kle, 'SSFEM');
            
            kle.getEffectiveKL();

            epsr = mc.MeshData.MatLibrary(mc.MeshData.TetType);

            %repeat epsr  Rn times, there are Rn columns, each column is a
            %sample
            epsr = repmat(epsr,[1 mc.Rn]);
            mc.SETUP.N = kle.getSdim();
            rng(mc.seed);
            mc.xi = randn(mc.Rn,length(kle.KLSet));

            for j = 1:length(kle.KLSet)
                KL = kle.KLSet{j};


                mat_idx = find(mc.MeshData.TetType == KL.domain);
                xi_i = mc.xi(:,j)'; % Rn columns, 1 Row
                xi_i = repmat(xi_i,[mc.MeshData.NT,1]); %Rn columns, NT rows
               
                if strcmp(level,'element')
                    K = KL.CORR.K_e;
                elseif strcmp(level,'full')
                    K = KL.CORR.K;
                end
                
                epsr(mat_idx,:) = epsr(mat_idx,:) + K.*xi_i(mat_idx,:)*KL.sd;
                
                mc.epsilon_r = epsr;
                mean_sd = mean(std(epsr(mat_idx,:)'));
                cprintf('*black','\tMean Standard Deviation of domain %d: %f\n',KL.domain, mean_sd);
            end

            mc.SETUP.domains = kle.getDomains();
            mc.SETUP.sds = kle.getSDs();
            mc.SETUP.nkls = kle.getNKLs();
            mc.SETUP.means = mc.getPermittivity(mc.SETUP.domains);
            mc.SETUP.type = 'MCS_EFF_KLE';
            mc.SETUP.p_order = -1;


        end
        
    
        function assignRandomMaterialVariation(mc, domains,sds)
            assignRandomMaterialVariation@SFEM(mc, domains,sds, 'MonteCarlo');
            rng(mc.seed);
            mc.xi = randn(mc.Rn,length(domains));
            epsr = mc.MeshData.MatLibrary(mc.MeshData.TetType);
            epsr = repmat(epsr,[1 mc.Rn]);
            for it=1:length(domains)
                domain = domains(it);
                sd = sds(it);
                mat_idx = find(mc.MeshData.TetType == domain);
                

                

                xi_i = repmat(mc.xi(:,it)',[length(mat_idx),1]);
                epsr(mat_idx,:) = epsr(mat_idx,:)+ sd*xi_i;
            end
            mc.epsilon_r = epsr;

            mc.SETUP.domains = domains;
            mc.SETUP.sds = sds;
            mc.SETUP.nkls = zeros(size(sds));
            mc.SETUP.means = mc.getPermittivity(mc.SETUP.domains);
            mc.SETUP.type = 'MCS_RV';

            mc.SETUP.p_order = -1;

        end


        function ret = assignSpatialMaterialVariation(mc, kle)
            %% Generate samples that are spatially correlated. The spatial
            %  correlation is represented by the KL Expansion


            assignSpatialMaterialVariation@SFEM(mc, kle, 'Monte Carlo');
            epsr = mc.MeshData.MatLibrary(mc.MeshData.TetType);

            %repeat epsr  Rn times, there are Rn columns, each column is a
            %sample
            epsr = repmat(epsr,[1 mc.Rn]);
            mc.SETUP.N = kle.getSdim();
            rng(mc.seed);
            mc.xi = randn(mc.Rn,mc.SETUP.N);

            for j = 1:length(kle.KLSet)
                KL = kle.KLSet{j};


                mat_idx = find(mc.MeshData.TetType == KL.domain);
               
                ret = zeros(100,4);
                %KL sum mean + Sigma_i=1^nkl sqrt(Lambda_i)*Phi_i*x_i;
                for i=1:KL.nkl
                    
                    k = (j-1)*KL.nkl + i;
%                     disp(k);
                    Phi_i = KL.CORR.Phi{i};  %1 column, NT rows
                    Phi_i = repmat(Phi_i,[1 mc.Rn]);  %Rn columns, NT rows
                    lambda_i = KL.CORR.Lambda(i);

                   
                    xi_i = mc.xi(:,k)'; % Rn columns, 1 Row
                    xi_i = repmat(xi_i,[mc.MeshData.NT,1]); %Rn columns, NT rows
                    tmp = epsr(mat_idx,:);
                    epsr = epsr + sqrt(lambda_i)*Phi_i.*xi_i*KL.sd;
                    ret(i,1) = max(mean(abs(tmp - 4.3)));
                    ret(i,2) = mean(mean(abs(tmp - 4.3)));
                    ret(i,3) = max(mean(abs(tmp - epsr(mat_idx,:))));
                    ret(i,4) = mean(mean(abs(tmp - epsr(mat_idx,:))));
                    ret(i,5) = vecnorm(mean(abs(tmp - 4.3)));
                    ret(i,6) =  mean(std((tmp - 4.3)));

                    fprintf('Max  diff for %d KLterms from det  = %f\n',i, max(mean(abs(tmp - 4.3))));
                    fprintf('Mean diff for %d KLterms from det  = %f\n',i, mean(mean(abs(tmp - 4.3))));
                    fprintf('Max  error for %d KLterms          = %f\n',i, max(mean(abs(tmp - epsr(mat_idx,:)))));
                    fprintf('Mean error for %d KLterms          = %f\n-----\n',i, mean(mean(abs(tmp - epsr(mat_idx,:)))));
                end
                mc.epsilon_r = epsr;
                mean_sd = mean(std(epsr(mat_idx,:)'));
                cprintf('*black','\tMean Standard Deviation of domain %d: %f\n',KL.domain, mean_sd);
            end

            mc.SETUP.domains = kle.getDomains();
            mc.SETUP.sds = kle.getSDs();
            mc.SETUP.nkls = kle.getNKLs();
            mc.SETUP.means = mc.getPermittivity(mc.SETUP.domains);
            mc.SETUP.type = 'MCS_KLE';
            mc.SETUP.p_order = -1;
        end
        
        function MCSimulation(mc,f)
            mc.calcElementalMatrix();
            mc.calcB();
            ref=zeros(mc.Rn,1);
            trans=zeros(mc.Rn,1);
            N = mc.Rn;
            WaitMessage = parfor_wait(N, 'Waitbar', true);
            parfor i=1:N
                WaitMessage.Send;
                mc.Assemble(mc.epsilon_r(:,i));
                mc.solve(f);
                ref(i)=mc.calcRef();
                trans(i)=mc.calcTrans();
            end
            WaitMessage.Destroy();
            mc.RESULTS.R=ref;
            mc.RESULTS.T=trans;
            mc.RESULTS.f = f;
            mc.RESULTS.seed = mc.seed;

           
            [a,b] = ksdensity(mc.RESULTS.T);
            mc.RESULTS.xpdf = b;
            mc.RESULTS.ypdf = a;
            mc.RESULTS.means21 = mean(mc.RESULTS.T);
            mc.RESULTS.stds21 = std(mc.RESULTS.T);
            
            mc.RESULTS.domains = mc.SETUP.domains ;
            mc.RESULTS.sds = mc.SETUP.sds;
            mc.RESULTS.means = mc.SETUP.means;
            mc.RESULTS.nkls = mc.SETUP.nkls;
            mc.RESULTS.p_order = mc.SETUP.p_order;

            mc.RESULTS.solvetime = -1;
            mc.RESULTS.type = mc.SETUP.type;
            mc.RESULTS.Rn = mc.Rn;
            mc.RESULTS.saved = false;
            mc.pushResult();
            mc.logresult();

        end
        
       
        
    end
end
