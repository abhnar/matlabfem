classdef StochasticCollocation< SFEM 
    properties
        
        
       
       
        epsilon_r;
        xi_sc;
        sc_base;
        scn;
        
       
        
       
        
    end
    

    methods
        
        function setSCN(obj,n)
            obj.scn =n;
        end
        function assignEffectiveMaterialVariation(sc, kle, level)
            if(nargin ==2) level = 'full';  end;

            assignEffectiveMaterialVariation@SFEM(sc, kle, 'SSFEM');
            
            kle.getEffectiveKL();

            epsr = sc.MeshData.MatLibrary(sc.MeshData.TetType);

            %repeat epsr  Rn times, there are Rn columns, each column is a
            %sample
            
            
            rng(sc.seed);
            sc.xi = randn(sc.Rn,length(kle.KLSet));
            N = length(kle.KLSet);
            sc.SETUP.N = N;
            rng(sc.seed);
            sc.xi = randn(sc.Rn,sc.SETUP.N);
            sc.sc_base = SC;
            sc.sc_base.setSDIM(N);
            sc.sc_base.setMeans(zeros(N,1));
            sc.sc_base.setSDs(ones(N,1));
            sc.sc_base.Rn = sc.Rn;
            sc.xi_sc  = sc.sc_base.generateCollocationPoints(sc.scn);
            sc.Rn = length(sc.xi_sc);
            epsr = repmat(epsr,[1 sc.Rn]);
            fprintf('Number of collocation points: %d\n',sc.Rn);

            for j = 1:length(kle.KLSet)
                KL = kle.KLSet{j};


                mat_idx = find(sc.MeshData.TetType == KL.domain);
                xi_i = sc.xi_sc(:,j)'; % Rn columns, 1 Row
                xi_i = repmat(xi_i,[sc.MeshData.NT,1]); %Rn columns, NT rows
               
                if strcmp(level,'element')
                    K = KL.CORR.K_e;
                elseif strcmp(level,'full')
                    K = KL.CORR.K;
                end
                
                epsr(mat_idx,:) = epsr(mat_idx,:) + K.*xi_i(mat_idx,:)*KL.sd;
                
                sc.epsilon_r = epsr;
                mean_sd = mean(std(epsr(mat_idx,:)'));
                cprintf('*black','\tMean Standard Deviation of domain %d: %f\n',KL.domain, mean_sd);
            end

            sc.SETUP.domains = kle.getDomains();
            sc.SETUP.sds = kle.getSDs();
            sc.SETUP.nkls = kle.getNKLs();
            sc.SETUP.means = sc.getPermittivity(sc.SETUP.domains);
            sc.SETUP.type = 'MCS_EFF_KLE';
            sc.SETUP.p_order = -1;


        end
        
    
        function assignRandomMaterialVariation(sc, domains,sds)
            assignRandomMaterialVariation@SFEM(sc, domains,sds, 'MonteCarlo');
            rng(sc.seed);
            sc.xi = randn(sc.Rn,length(domains));
            epsr = sc.MeshData.MatLibrary(sc.MeshData.TetType);
            
            sc.SETUP.N = length(domains);
            rng(sc.seed);
            sc.xi = randn(sc.Rn,sc.SETUP.N);
            sc.sc_base = SC;
            sc.sc_base.setSDIM(sc.SETUP.N);
            sc.sc_base.setMeans(zeros(sc.SETUP.N,1));
            sc.sc_base.setSDs(ones(sc.SETUP.N,1));
            sc.sc_base.Rn = sc.Rn;
            sc.xi_sc  = sc.sc_base.generateCollocationPoints(sc.scn);
            sc.Rn = length(sc.xi_sc);
            epsr = repmat(epsr,[1 sc.Rn]);
            fprintf('Number of collocation points: %d\n',sc.Rn);

            epsr = repmat(epsr,[1 sc.Rn]);
            for it=1:length(domains)
                domain = domains(it);
                sd = sds(it);
                mat_idx = find(sc.MeshData.TetType == domain);
                

                

                xi_i = repmat(sc.xi_sc(:,it)',[length(mat_idx),1]);
                epsr(mat_idx,:) = epsr(mat_idx,:)+ sd*xi_i;
            end
            sc.epsilon_r = epsr;

            sc.SETUP.domains = domains;
            sc.SETUP.sds = sds;
            sc.SETUP.nkls = zeros(size(sds));
            sc.SETUP.means = sc.getPermittivity(sc.SETUP.domains);
            sc.SETUP.type = 'MCS_RV';

            sc.SETUP.p_order = -1;

        end


        function ret = assignSpatialMaterialVariation(sc, kle)
            %% Generate samples that are spatially correlated. The spatial
            %  correlation is represented by the KL Expansion


            assignSpatialMaterialVariation@SFEM(sc, kle, 'Stochastic Collocation');
            epsr = sc.MeshData.MatLibrary(sc.MeshData.TetType);

            %repeat epsr  Rn times, there are Rn columns, each column is a
            %sample
            
            sc.SETUP.N = kle.getSdim();
            rng(sc.seed);
            sc.xi = randn(sc.Rn,sc.SETUP.N);
            sc.sc_base = SC;
            sc.sc_base.setSDIM(sc.SETUP.N);
            sc.sc_base.setMeans(zeros(sc.SETUP.N,1));
            sc.sc_base.setSDs(ones(sc.SETUP.N,1));
            sc.sc_base.Rn = sc.Rn;
            sc.xi_sc  = sc.sc_base.generateCollocationPoints(sc.scn);
            sc.Rn = length(sc.xi_sc);
            epsr = repmat(epsr,[1 sc.Rn]);
            

            for j = 1:length(kle.KLSet)
                KL = kle.KLSet{j};


                mat_idx = find(sc.MeshData.TetType == KL.domain);
               
                ret = zeros(100,4);
                %KL sum mean + Sigma_i=1^nkl sqrt(Lambda_i)*Phi_i*x_i;
                for i=1:KL.nkl
                    
                    k = (j-1)*KL.nkl + i;
%                     disp(k);
                    Phi_i = KL.CORR.Phi{i};  %1 column, NT rows
                    Phi_i = repmat(Phi_i,[1 sc.Rn]);  %Rn columns, NT rows
                    lambda_i = KL.CORR.Lambda(i);

                   
                    xi_i = sc.xi_sc(:,k)'; % Rn columns, 1 Row
                    xi_i = repmat(xi_i,[sc.MeshData.NT,1]); %Rn columns, NT rows
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
                sc.epsilon_r = epsr;
                mean_sd = mean(std(epsr(mat_idx,:)'));
                cprintf('*black','\tMean Standard Deviation of domain %d: %f\n',KL.domain, mean_sd);
            end

            sc.SETUP.domains = kle.getDomains();
            sc.SETUP.sds = kle.getSDs();
            sc.SETUP.nkls = kle.getNKLs();
            sc.SETUP.means = sc.getPermittivity(sc.SETUP.domains);
            sc.SETUP.type = 'MCS_KLE';
            sc.SETUP.p_order = -1;
        end
        
        function SCSimulation(sc,f)
            sc.calcElementalMatrix();
            sc.calcB();
            ref=zeros(sc.Rn,1);
            trans=zeros(sc.Rn,1);
            N = sc.Rn;
            WaitMessage = parfor_wait(N, 'Waitbar', true);
            parfor i=1:N
                WaitMessage.Send;
                sc.Assemble(sc.epsilon_r(:,i));
                sc.solve(f);
                ref(i)=sc.calcRef();
                trans(i)=sc.calcTrans();
            end
            WaitMessage.Destroy();
            sc.RESULTS.R=ref;
            sc.RESULTS.T=trans;
            sc.RESULTS.T = sc.sc_base.collocate(sc.RESULTS.T);
            sc.RESULTS.R = sc.sc_base.collocate(sc.RESULTS.R);
            sc.RESULTS.f = f;
            sc.RESULTS.seed = sc.seed;

           
            [a,b] = ksdensity(sc.RESULTS.T);
            sc.RESULTS.xpdf = b;
            sc.RESULTS.ypdf = a;
            sc.RESULTS.means21 = mean(sc.RESULTS.T);
            sc.RESULTS.stds21 = std(sc.RESULTS.T);
            
            sc.RESULTS.domains = sc.SETUP.domains ;
            sc.RESULTS.sds = sc.SETUP.sds;
            sc.RESULTS.means = sc.SETUP.means;
            sc.RESULTS.nkls = sc.SETUP.nkls;
            sc.RESULTS.p_order = sc.SETUP.p_order;

            sc.RESULTS.solvetime = -1;
            sc.RESULTS.type = sc.SETUP.type;
            sc.RESULTS.Rn = sc.Rn;
            sc.RESULTS.saved = false;
            sc.pushResult();
            sc.logresult();

        end
        
       
        
    end
end
