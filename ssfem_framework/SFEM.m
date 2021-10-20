classdef SFEM < handle & EFEM
    properties
        Rn;
        iskl = false;
        RESULTS;
        CACHE;
        seed = 1;
        RV;
        KLE;
        nkl = -1;
        p_order = -1;
        xi;
        rescnt = 0;
        
    end
    
     properties(Access = protected)
        linestyle = {'-','--''-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':'};
        linewidth = 1.5;
        markers = {'o','o','o','o','s','s','s','s','*','*','*','*','v','v','v','v'};
        colors = {'blue','blue','blue','blue','black','black','black','black','red','red','red','red','green','green','green','green'};
        
     end
    
    methods
        
        function init(mc, R, varargin)
            %% R = number of random samples
            % 'kle', true <- enables KLE
            % 'nklterms', numberofklterms
            %
            mc.iskl = false;
            mc.Rn = R;
            mc.p_order = -1;
            mc.nkl = -1;
            if(length(varargin) > 1 )
                for i=1:2:length(varargin)
                    switch(varargin{i})
                        case 'kle'
                            
                            mc.iskl = varargin{i+1};
                        case 'nklterms'
                            if(~isnumeric(varargin{i+1}))
                                error('Expection a numeric value for nklterms');
                            end
                            mc.nkl = varargin{i+1};
                        case 'p_order'
                            
                            mc.p_order = varargin{i+1};
                            
                            
                    end
                end
                
            end
            if(mc.nkl == -1 && mc.iskl)
                fprintf("Setting default value of nkl = 3\n");
                mc.nkl = 3;
            end
            if(mc.p_order == -1)
                fprintf("Setting default value of p_order = 2\n");
                mc.p_order = 2;
            end
            
        end
        
        
        function assignRandomMaterialVariation(ssfem, domain, sd, type)
            cprintf('*blue','----------------------------------------\n');
            str = sprintf('         %s (RV - Fixed over Space)         \n', type);
            cprintf('*blue',str);
            cprintf('*blue','----------------------------------------\n');
            cprintf('*black','\tRandom seed             : %d\n', ssfem.seed);
            cprintf('*black','\tDomain                  : %d\n', domain);
            cprintf('*black','\tInput Standard Deviation: %2.2f\n', sd);
        end
        
        
        function assignSpatialMaterialVariation(ssfem, kle, type)
            %% Generate samples that are spatially correlated. The spatial
            %  correlation is represented by the KL Expansion
            kldata = kle.CORR;
            cprintf('*blue','-----------------------------------------------------\n');
            cprintf('*blue','         %s Simulation (KLE - Spatial)         \n', type);
            cprintf('*blue','-----------------------------------------------------\n');
          
            
            %cprintf('_blue','                                         \n');
            cprintf('*black','\tRandom seed       : %d\n', ssfem.seed);
            cprintf('_blue','KL Information\n');
            cprintf('*black','\tDomain            : %d\n', kldata.domain);
            cprintf('*black','\tCorrelation Length: %e m (factor = %2.2f)\n', kldata.corlen, kle.KLDATA.corlenfactor);
            cprintf('*black','\tNumber of KL Terms: %d\n', ssfem.nkl);
        end
        
        function plot(ssfem,type)
            if nargin <2
                fprintf('Argument should be a string in ("fsweep","current")\n')
                return
            end
            if(strcmpi(type,'fsweep'))
                results = ssfem.CACHE;
                if length(results) < 2
                    fprintf('At least two results should be in cache');
                    return;
                end
                plot_samples = zeros(ssfem.Rn,length(results));
                freq = [];
                for i=1:length(results)
                    plot_samples(:,i) = results{i}.T;
                    freq(i) = results{i}.f;
                end
                
                for j = 1:100
                    plot(freq,plot_samples(j,:));
                    hold on;
                end
                xlabel('Frequency (GHz)')
                ylabel('|S_{21}|');
            elseif (strcmpi(type,'current'))
                [a,b] = ksdensity(ssfem.RESULTS.T);
                plot(b,a)
            end
        end
        function pushResult(mc)
            
            mc.CACHE{mc.rescnt+1} =  mc.RESULTS;
            mc.rescnt = mc.rescnt +1;
        end
        
        function clearCache(mc)
            mc.rescnt = 0;
            mc.CACHE = {};
        end
        
        function saveCache(mc, fname)
            CACHE = mc.CACHE;
            save([fname,'.mat'],'CACHE');
        end
        function loadCache(mc, fname)
            
            load([fname,'.mat']);
            mc.CACHE = CACHE;
        end
        function setSeed(mc, sd)
            if(nargin == 1)
                mc.seed = round(randn*100,0);
            else
                mc.seed =  sd;
            end
            rng(mc.seed);
        end
        
        function val = getNorm(mc,kle,i1)
            
            supSet = kle.KLDATA.supSet;
            a=supSet.Co(:,:,1);
            b=supSet.Co(:,:,2);
            c=supSet.Co(:,:,3);
            d=supSet.Co(:,:,4);
            V=supSet.V;
            
            
            tet=supSet.tet;
            
            
            Evec = kle.CORR.Evec;
            
            ix1 = supSet.ix1;
            nds = supSet.nds;
            
            
            
            eigvec1=Evec(:,i1);
            eigvec2=Evec(:,i1);
            val=0;
            for i=1:length(tet)
                
                %tetid=r(i);
                
                lix = i;
                Xc=mean(nds(tet(lix,:),:));
                
                
                integ=0;
                for k=1:4
                    Lk = @(x) 1/(6*V(lix))*(a(lix,k)+b(lix,k)*x(1) +c(lix,k)*x(2)+d(lix,k)*x(3));
                    
                    integ = integ+Lk(Xc)*eigvec1(tet(lix,k))*eigvec2(tet(lix,k));
                end
                
                
                val=val+integ*V(lix);
            end
            
            
            
            val = sqrt(val);
            
        end
    end
end