classdef SFEM < handle & EFEM
    properties
        Rn;
        iskl = false;
        RESULTS;
        CACHE;
        seed = 1;
        RV;
        SETUP;
        nkl = -1;
        p_order = -1;
        xi;
        rescnt = 0;
        ID = -1;
        
    end
    
     properties(Access = protected)
        linestyle = {'-','--''-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':'};
        linewidth = 1.5;
        markers = {'o','o','o','o','s','s','s','s','*','*','*','*','v','v','v','v'};
        colors = {'blue','blue','blue','blue','black','black','black','black','red','red','red','red','green','green','green','green'};
        
     end
    
    methods
        function mats = getPermittivity(sfem, doms)
            mats = sfem.MeshData.MatLibrary(doms);
        end
        function clear(sfem)
            sfem.CAHCE = [];
            sfem.d_stoch = [];
            sfem.RESULTS = [];
            sfem.SETUP = [];
            sfem.B= [];
            sfem.RHS= [];
            sfem.E_freespace= [];
            sfem.F_freespace= [];
            sfem.E_material= [];
            sfem.F_material= [];
            sfem.EIGSYS = [];
            sfem.temp = [];
            sfem.ipfacefield = [];
            sfem.pecedge = [];
            sfem.EM = [];
            sfem.MeshData = [];
            sfem.xi = [];
        end
        function ssfem = SFEM()
            x = 0;while(x<10000)x = round(rand*100000); end; ssfem.ID = x;
             
        end
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
            
            cprintf('*Blue','-----------------------------------------------------\n');
            cprintf('*Blue','         %s Simulation (KLE - Spatial)         \n', type);
            cprintf('*Blue','-----------------------------------------------------\n');
            
            ndom = length(kle.KLSet);
            cprintf('*black','\tRandom seed       : %d\n', ssfem.seed);
            cprintf('*black','\tNumber of Domains : %d\n', ndom);
            for i=1:ndom
                cprintf('*black','\t\tDomain : %d\n', kle.KLSet{i}.domain);
                cprintf('*black','\t\t\tStandard Deviation     : %2.2f\n', kle.KLSet{i}.sd);
                cprintf('*black','\t\t\tCorrelation Length     : %e m (factor = %2.2f)\n', kle.KLSet{i}.CORR.corlen, kle.KLSet{i}.corlenfactor);
                cprintf('*black','\t\t\tNumber of KL Terms     : %d\n', kle.KLSet{i}.nkl);
            end
            
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
        function logresult(sfem)
            if ~exist('.results', 'dir')
                mkdir('.results');
            end
            timestamp = datestr(datetime);
            doms = strtrim(sprintf('%d ',sfem.RESULTS.domains));
            sds = strtrim(sprintf('%f ',sfem.RESULTS.sds));
            nkls = strtrim(sprintf('%d ',sfem.RESULTS.nkls));
            means = strtrim(sprintf('%f ', real(sfem.RESULTS.means)));
            res= sprintf('%s %f %d %f %f %d (%s) (%s) (%s) (%s) %f', sfem.RESULTS.type, sfem.RESULTS.f,sfem.RESULTS.Rn, ...
                sfem.RESULTS.means21,sfem.RESULTS.stds21,sfem.RESULTS.p_order,doms, means, sds, nkls,sfem.RESULTS.solvetime);
            id = fopen('.results/result_log.txt','a+');
%             disp(id)
            fwrite(id, sprintf('%s %d  %s\n',timestamp, sfem.ID, res));
            fclose(id);
        end
        function pushResult(mc)
            fname = ['.results/',strrep(strrep(datestr(datetime),' ','_'),':','-'),'.mat'];
            
            if ~exist('.results', 'dir')
                mkdir('.results');
            end
            
            if ~exist('.results/result_cnt.txt','file')
%                 disp('Fresh Result Set');
                id = fopen('.results/result_cnt.txt','w');
                fwrite(id,'0');
                fclose(id);
                
            end

            id = fopen('.results/result_cnt.txt','r');
            x = str2double(fgetl(id));
%             disp(x);
            fclose(id);
           
            id = fopen('.results/result_table.txt','a+');
%             disp(id)
            fwrite(id, sprintf('%d %d %s N\n',mc.ID, x+1,fname));
            fclose(id);

            id = fopen('.results/result_cnt.txt','w');
            fwrite(id,num2str(x+1) );
            fclose(id);
            mc.RESULTS.saved = true;
            RESULTS = mc.RESULTS;
            save(fname,'RESULTS','-mat');
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
     
    end
end