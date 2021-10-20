classdef MonteCarlo< SFEM 
    properties
        
        
       
       
        epsilon_r;
        
       
        
       
        
    end
    

    methods
        
        
       
        
    
        function assignRandomMaterialVariation(mc, domain,sd)
            assignRandomMaterialVariation@SFEM(mc, domain,sd, 'MonteCarlo');
           
            mat_idx = find(mc.MeshData.TetType == domain);
            epsr = mc.MeshData.MatLibrary(mc.MeshData.TetType);
            epsr = repmat(epsr,[1 mc.Rn]);
            
            rng(mc.seed);
            mc.xi = randn(mc.Rn,1);
            xi_i = repmat(mc.xi',[length(mat_idx),1]);
            epsr(mat_idx,:) = epsr(mat_idx,:)+ sd*xi_i;
            mc.epsilon_r = epsr;
            
        end
        
        
        function assignSpatialMaterialVariation(mc, kle)
        %% Generate samples that are spatially correlated. The spatial 
        %  correlation is represented by the KL Expansion
        
        
        assignSpatialMaterialVariation@SFEM(mc, kle, 'Monte Carlo');
        
        kldata = kle.CORR;
        
        mat_idx = find(mc.MeshData.TetType == kldata.domain);
        epsr = mc.MeshData.MatLibrary(mc.MeshData.TetType);
        
        %repeat epsr  Rn times, there are Rn columns, each column is a
        %sample
        epsr = repmat(epsr,[1 mc.Rn]);
        rng(mc.seed);
        mc.xi = randn(mc.Rn,mc.nkl);
        %KL sum mean + Sigma_i=1^nkl sqrt(Lambda_i)*Phi_i*x_i;
        for i=1:mc.nkl
            Phi_i = kldata.Phi{i};  %1 column, NT rows
            Phi_i = repmat(Phi_i,[1 mc.Rn]);  %Rn columns, NT rows
            lambda_i = kldata.Lambda(i);
            norm_i = kldata.Norms(i);
            norm_i = mc.getNorm(kle,i);
            xi_i = mc.xi(:,i)'; % Rn columns, 1 Row
            xi_i = repmat(xi_i,[mc.MeshData.NT,1]); %Rn columns, NT rows
            epsr = epsr + sqrt(lambda_i)*Phi_i.*xi_i*kle.KLDATA.sd/norm_i;
        end
        mc.epsilon_r = epsr;
        mean_sd = mean(std(epsr(mat_idx,:)'));
        cprintf('*black','\tMean Standard Deviation: %f\n', mean_sd);
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
        
        end
        
       
        
    end
end
