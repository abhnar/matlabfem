classdef SSFEM<  SFEM
    properties

        d_stoch;
        show_progress = true;
        loadedresult;
        
    end
   
    methods
        function analyseresult(ssfem,cmd, id)
            if strcmp(cmd,'load')
                id = fopen('results/femlog.txt');
                line = fgetl(id);
                ids = [];
                res = {};
                cnt = 1;
                idmap = [];
                fmap = {};
                type = {};
                while(line~=-1)
                    tmp = strsplit(line, ' ');
                    if ~(ismember(str2double(tmp{1}), ids))
                        ids = [ids str2double(tmp{1})];
                    end
                    
                    res{cnt} = line;
                    idmap(cnt) = str2double(tmp{1});
                    t = split(tmp{4},'/');t = t{end};
                    type{cnt} = tmp{5};
                    fmap{cnt} = replace(['results/', t,tmp{2},'_', tmp{3},'.mat'],':','-');
                    cnt = cnt + 1;
                    line = fgetl(id);
                    
                    
                end
                fclose(id);
                ssfem.loadedresult.ids = ids;
                ssfem.loadedresult.res = res;
                ssfem.loadedresult.fmap = fmap;
                ssfem.loadedresult.idmap = idmap;
            end
            if strcmp(cmd,'list')
                disp(ssfem.loadedresult.ids');
            end
            if strcmp(cmd,'band')
                cnt = 0;
                flist = ssfem.loadedresult.fmap(ssfem.loadedresult.idmap == id);
                for i=1:length(flist)
                    load(flist{i});
                    [y,x] = ksdensity(T);
                    cnt = cnt+1;
                    mn(cnt) = mean(T);
                    Tu = T(T>=mn(cnt));
                    Td = T(T<=mn(cnt));
                    sd(cnt) = std(T);
                    sdu(cnt) = std(Tu);
                    sdd(cnt) = std(Td);
                    freq(cnt) = f;
                end
                [freq,ix] = sort(freq);
                mn = mn(ix);
                sdu = sdu(ix);
                sdd=sdd(ix);
                hold on;
                patch([freq fliplr(freq)], [mn+sdu*3 fliplr(mn-sdd*3)], 'g', 'facealpha', 0.5)
                plot(freq,mn);
                ssfem.loadedresult.freq = freq;
                
            end

             if strcmp(cmd,'table')
                cnt = 0;
                flist = ssfem.loadedresult.fmap(ssfem.loadedresult.idmap == id);
                tps = cell(length(flist),1);
                for i=1:length(flist)
                    load(flist{i});
                    [y,x] = ksdensity(T);
                    cnt = cnt+1;
                    mn(cnt) = mean(T);
                    Tu = T(T>=mn(cnt));
                    Td = T(T<=mn(cnt));
                    sd(cnt) = std(T);
                    sdu(cnt) = std(Tu);
                    sdd(cnt) = std(Td);
                    freq(cnt) = f;
                    tps{cnt} = type;
                end
                [freq,ix] = sort(freq);
                mn = mn(ix);
                sdu = sdu(ix);
                sdd=sdd(ix);
                sd=sd(ix);
                fprintf("slno\tfrequency\t   Mean\t     std\n")
                fprintf('-----------------------------------------\n')
                for i =1:length(freq)
                    fprintf("%4d\t%3.6f\t%1.6f\t%1.6f\t%s\n", i, freq(i),mn(i),sd(i), tps{i});
                end
                
            end
        
        end
        
        function execute(ssfem,fn, id)
        if(nargin ==3)
          exec_id = id;
        else
          exec_id = round(rand*1e5);  
        end
        fprintf('Execution id: %d\n', exec_id);
            
            id = fopen(fn);
            line = fgetl(id);
            status = 0;
            while(line~=-1)
%                 disp(line);
                if status == 0
                    if(strcmp(line,'BEGIN'))
                        status = 1;
                    end
                else
                    model = split(line,':');
                    model = model{2};
                    
                    mats = split(fgetl(id),':');
                    mats = split(strtrim(mats{2}),',');
                    mats = str2double(mats);
                    
                    doms = split(fgetl(id),':');
                    doms = split(strtrim(doms{2}),',');
                    doms = str2double(doms);
                    
                    stds = split(fgetl(id),':');
                    stds = split(strtrim(stds{2}),',');
                    stds = str2double(stds);
                    
                    klns = split(fgetl(id),':');
                    klns = split(strtrim(klns{2}),',');
                    klns = str2double(klns);
                    
                    corlen = split(fgetl(id),':');
                    corlen = str2double(strtrim(corlen{2}));
                    
                    p = split(fgetl(id),':');
                    p = str2double(strtrim(p{2}));
                    
                    freq = split(fgetl(id),':');
                    
                    if(strcmp('r',freq{2}))
                        freq = split(freq{3},',');
                        freq = linspace(str2double(freq{1}),str2double(freq{2}),str2double(freq{3}));
                    elseif(strcmp('d',freq{2}))
                         freq = split(freq{3},',');
                         freq = str2double(freq);
                    else
                        error("Freq type should be either 'r' or 'd'");
                    end
                   
                    
                    seed = split(fgetl(id),':');
                    seed = str2double(strtrim(seed{2}));
                    
                    Rn = split(fgetl(id),':');
                    Rn = str2double(strtrim(Rn{2}));
                    
                    type = split(fgetl(id),':');
                    type = strtrim(type{2});
                    
                    tag = split(fgetl(id),':');
                    tag = strtrim(tag{2});
                    status = 0;
                    
                    ssfem.verbose = 0;
                    ssfem.LoadMesh(model);
                    ssfem.meshStatistics();
                   
                    ssfem.setMaterials(mats);
                    ssfem.buildSystem();
                    
                    
                    
                    %% KLE Setup
                    kle = KLE3D;
                    kle.process(ssfem,doms);
                    kle.evaluate_Phi('nlambda',max(klns),'corlens',corlen);
                    kle.getKLEData(corlen);
                    kle.setSD(stds);
                    kle.setNKL(klns);
                    ssfem.init(Rn,'kle',true,'p_order',p);
                    
                    ssfem.setSeed(seed);
                    if strcmp(type,'kle')
                        ssfem.assignSpatialMaterialVariation(kle);
                    elseif strcmp(type,'kle_eff')
                        kle.getEffectiveKL();
                        ssfem.assignEffectiveMaterialVariation(kle);
                    end
                    for i =1:length(freq)
                        f = freq(i);
                       
                        ssfem.ssfemsimulation(f);
                        dt = datestr(datetime);
                        s = sprintf('%5d %20s %20s %10s %15s %15s %15s %2.2f %10s %2d %3d %3.6f %2.6f %2.6f %3.6f %s\n',exec_id,dt, ssfem.Meshfile,type,['[',sprintf('%d,',doms),']'],...
                           ['[',sprintf('%2.2f,',real(ssfem.getPermittivity(doms))),']'], ['[',sprintf('%2.3f,',stds),']'], corlen, ['[',sprintf('%d,',klns),']'], ...
                            p, ssfem.RESULTS.P,f,ssfem.RESULTS.means21, ssfem.RESULTS.stds21,  ssfem.RESULTS.solvetime, tag);
                        logfid = fopen('results/femlog.txt','a+');
                        fwrite(logfid,s);
                        fclose(logfid);
                        
                        fname = split(ssfem.Meshfile,'/');
                        fname = ['results/', fname{end}, strrep(strrep(dt,':','-'),' ','_') ,'.mat'];
                        T = ssfem.RESULTS.T;
                        [pdfy,pdfx] = ksdensity(ssfem.RESULTS.T);
                        save(fname, 'pdfx','pdfy','T','f','dt','p','type','exec_id','tag');
                        ssfem.clearresults();

                    end
                    
                    
                    
                end
                line = fgetl(id);
            end
            
            fclose(id);
        
        end
        function clearresults(ssfem)
            ssfem.CACHE = [];
            ssfem.RESULTS = [];
            ssfem.SETUP.Et = [];
        end
        function assignEffectiveMaterialVariation(ssfem, kle, level)
            if(nargin ==2) 
                level = 'full'; 
            end

            assignEffectiveMaterialVariation@SFEM(ssfem, kle, 'SSFEM');
            domains = kle.getDomains();
            
            sds = kle.getSDs();
            kle.getEffectiveKL();
            sds_ = kle.getEffSD();

            ssfem.SETUP.TSet = cell(length(domains),1);
            for it = 1:length(domains)
                domain = domains(it);
                

                ssfem.buildSystem();
                
               
                
                temp=ssfem.MeshData.TetType;
                temp(temp~=domain)=0;
                temp(temp==domain)=1;
                if(strcmp(level,'full'))
                    sd  = sds_(it);
                    epr = sd.*(temp);
                elseif(strcmp(level,'element'))
                    sd = sds(it);
                    epr = temp;
                    epr = epr.*kle.KLSet{it}.CORR.K_e*sd;
                end
                
                ssfem.Assemble(epr);
                
                ssfem.SETUP.Tset{it+1} = ssfem.T;
                
                ssfem.SETUP.Tset{it+1}(ssfem.pecedge,:) = [];
                ssfem.SETUP.Tset{it+1}(:,ssfem.pecedge) = [];
                
            end
            ssfem.SETUP.N = length(domains);
            ssfem.SETUP.P = PC.getP(ssfem.SETUP.N,ssfem.p_order);
            ssfem.SETUP.cijk = PC.c_ijk(ssfem.SETUP.N,ssfem.p_order,1);
            ssfem.SETUP.Psi = PC.getHermite_PC(ssfem.SETUP.N,ssfem.p_order);


            ssfem.SETUP.domains = domains;
            ssfem.SETUP.means = ssfem.getPermittivity(ssfem.SETUP.domains);
            ssfem.SETUP.sds = sds;
            ssfem.SETUP.type = 'SSFEM_EFF_KLE';

            ssfem.SETUP.nkls = zeros(size(sds));
            ssfem.SETUP.p_order = ssfem.p_order;

            rng(ssfem.seed);
            ssfem.xi = randn(ssfem.Rn,ssfem.SETUP.N);

            fprintf('\t');
            cprintf('_black','Stochastic Setup\n');
            cprintf('*black','\t\tN : %d\n', ssfem.SETUP.N);
            cprintf('*black','\t\tP : %d\n', ssfem.SETUP.P);
 
        end
      
        
        function assignRandomMaterialVariation(ssfem, domains, sds)
            
            assignRandomMaterialVariation@SFEM(ssfem, domains, sds, 'SSFEM');
            ssfem.SETUP.TSet = cell(length(domains),1);
            for it = 1:length(domains)
                domain = domains(it);
                sd = sds(it);

                ssfem.buildSystem();
                
               
                
                temp=ssfem.MeshData.TetType;
                temp(temp~=domain)=0;
                temp(temp==domain)=1;
                epr = sd.*(temp);
                
                ssfem.Assemble(epr);
                
                ssfem.SETUP.Tset{it+1} = ssfem.T;
                
                ssfem.SETUP.Tset{it+1}(ssfem.pecedge,:) = [];
                ssfem.SETUP.Tset{it+1}(:,ssfem.pecedge) = [];
                
            end
            ssfem.SETUP.N = length(domains);
            ssfem.SETUP.P = PC.getP(ssfem.SETUP.N,ssfem.p_order);
            ssfem.SETUP.cijk = PC.c_ijk(ssfem.SETUP.N,ssfem.p_order,1);
            ssfem.SETUP.Psi = PC.getHermite_PC(ssfem.SETUP.N,ssfem.p_order);


            ssfem.SETUP.domains = domains;
            ssfem.SETUP.means = ssfem.getPermittivity(ssfem.SETUP.domains);
            ssfem.SETUP.sds = sds;
            ssfem.SETUP.type = 'SSFEM_RV';

            ssfem.SETUP.nkls = zeros(size(sds));
            ssfem.SETUP.p_order = ssfem.p_order;

            rng(ssfem.seed);
            ssfem.xi = randn(ssfem.Rn,ssfem.SETUP.N);

            fprintf('\t');
            cprintf('_black','Stochastic Setup\n');
            cprintf('*black','\t\tN : %d\n', ssfem.SETUP.N);
            cprintf('*black','\t\tP : %d\n', ssfem.SETUP.P);
 
        end
        
        function assignSpatialMaterialVariation(ssfem, kle)
            
             assignSpatialMaterialVariation@SFEM(ssfem, kle, 'SSFEM');
            
            
            
            ssfem.buildSystem();
            
            
            for j=1:length(kle.KLSet)
                KL = kle.KLSet{j};
                for i=1:KL.nkl
                    k = (j-1)*KL.nkl + i +1;
%                     disp(k)
                    
                    epr = KL.sd*sqrt(KL.CORR.Lambda(i))*KL.CORR.Phi{i};
                    ssfem.Assemble(epr);
                    ssfem.SETUP.Tset{k} =  ssfem.T;
                    ssfem.SETUP.Tset{k}(:,ssfem.pecedge)=[];
                    ssfem.SETUP.Tset{k}(ssfem.pecedge,:)=[];
                end
            end

            ssfem.SETUP.N = kle.getSdim();
            ssfem.SETUP.P = PC.getP(ssfem.SETUP.N,ssfem.p_order);
            ssfem.SETUP.cijk = PC.c_ijk(ssfem.SETUP.N,ssfem.p_order,1);
            ssfem.SETUP.Psi = PC.getHermite_PC(ssfem.SETUP.N,ssfem.p_order);

            ssfem.SETUP.domains = kle.getDomains();
            ssfem.SETUP.means = ssfem.getPermittivity(ssfem.SETUP.domains);
            ssfem.SETUP.sds = kle.getSDs();
            ssfem.SETUP.nkls = kle.getNKLs();
            ssfem.SETUP.type = 'SSFEM_KLE';
            ssfem.SETUP.p_order = ssfem.p_order;

            rng(ssfem.seed);
            ssfem.xi = randn(ssfem.Rn,ssfem.SETUP.N);

            fprintf('\t');
            cprintf('_black','Stochastic Setup\n');
            cprintf('*black','\t\tN : %d\n', ssfem.SETUP.N);
            cprintf('*black','\t\tP : %d\n', ssfem.SETUP.P);
        
        end
        
        function ssfemsimulation(ssfem,f)
            ssfem.buildSystem();
            [K, b]=ssfem.buildDeterministicSystem(f);
            
            Tset = ssfem.SETUP.Tset;
            Tset{1} = K;
            for i=2:length(Tset)
                Tset{i} = -ssfem.K0^2*Tset{i};
            end

            cijk = ssfem.SETUP.cijk;
            P = ssfem.SETUP.P;
            
            
            
            Ndof = size(K,1);
            ssfem.SETUP.Ndof = Ndof;
            X = cell(ssfem.SETUP.P,ssfem.SETUP.P);
            Y = cell(ssfem.SETUP.P,ssfem.SETUP.P);
            VALS = cell(ssfem.SETUP.P,ssfem.SETUP.P);

            TSet = Tset;
            N = ssfem.SETUP.N;
            P = ssfem.SETUP.P;
           
            
            
           
            G=zeros(Ndof*P,1);
            G(1:Ndof)=b;
            
%            save('C:\Users\abhijith\workspace\python\python_ssfem_with_MATLAB\SETUP.mat','TSet','G','cijk','P','N')

            if(ssfem.show_progress)
                wb = waitbar(0,"Assembling SSFEM",'Name', 'SSFEM');
            end

            for j = 1 : P
                
                for k = j : P
                    S = sparse(Ndof,Ndof);
                    
                    for i = 1 : ssfem.SETUP.N + 1
                        S =  S + cijk{i}(j,k) * Tset{i};
                    end
                    
                    [x,y] = find(S);
                    val = S(S~=0);
                    X{j,k} = x;
                    Y{j,k} = y;
                    VALS{j,k} = val;
                    
                    
                end
                
                if(ssfem.show_progress)
                    waitbar(j/P,wb,"Assembling SSFEM")
                end
                
            end
            
            for j = 1 : P
                for k = 1 : j
                    X{j,k} = X{k,j} ;
                    Y{j,k} = Y{k,j} ;
                    VALS{j,k} = VALS{k,j};
                end
            end
            for j = 1 : P
                for k = 1 : P
                    X{j,k} = X{j,k} + (j-1)*Ndof;
                    Y{j,k} = Y{j,k} + (k-1)*Ndof;
                    
                end
            end

            if(ssfem.show_progress)
                close(wb);
            end

            X = vertcat(X{:});
            Y = vertcat(Y{:});
            VALS= vertcat(VALS{:});
            ST = sparse(X,Y,VALS);
            ssfem.RESULTS.NZ = (nnz(ST));
            
            clear X Y VALS S x y val;
            
            
            fprintf('Solving at frequency %f\n',f);
            tic;
            xr = ST\G;
            t1 = toc;
          
            ssfem.d_stoch = reshape(xr,Ndof,P);
            
            
            ssfem.d_stoch = reshape(xr,Ndof,P);
            
            
            ssfem.EvaluatePCE();

            ssfem.RESULTS.T = ssfem.calcTrans_Stoch(ssfem.SETUP.Et);
            [a,b] = ksdensity(ssfem.RESULTS.T);
            ssfem.RESULTS.xpdf = b;
            ssfem.RESULTS.ypdf = a;
            ssfem.RESULTS.means21 = mean(ssfem.RESULTS.T);
            ssfem.RESULTS.stds21 = std(ssfem.RESULTS.T);
            
            ssfem.RESULTS.domains = ssfem.SETUP.domains ;
            ssfem.RESULTS.sds = ssfem.SETUP.sds;
            ssfem.RESULTS.nkls = ssfem.SETUP.nkls;
            ssfem.RESULTS.solvetime = t1;
            disp(t1);
%             ssfem.RESULTS.N = size(ST,1);
%             ssfem.RESULTS.NZ = nnz(ST);
            ssfem.RESULTS.P=P;
            ssfem.RESULTS.p_order=ssfem.SETUP.p_order;
            
            ssfem.RESULTS.means = ssfem.SETUP.means;
            ssfem.RESULTS.type = ssfem.SETUP.type;
            ssfem.RESULTS.f = f;
            ssfem.RESULTS.seed = ssfem.seed;
            ssfem.RESULTS.Rn = ssfem.Rn;
            ssfem.RESULTS.saved = false;
            ssfem.pushResult();
            ssfem.logresult();
        end

        function resample(ssfem, seed, Rn)
            ssfem.RESULTS.seed = seed;
            ssfem.Rn = Rn;
            rng(seed);
            ssfem.xi = randn(Rn, ssfem.SETUP.N);
            ssfem.EvaluatePCE();
            ssfem.RESULTS.T = ssfem.calcTrans_Stoch(ssfem.SETUP.Et);
            
            ssfem.RESULTS.seed = ssfem.seed;
            ssfem.pushResult();
            ssfem.logresult();
        end

        function EvaluatePCE(ssfem)
            disp('Evaluating PCE');
            d = ssfem.d_stoch;

            Psi=ssfem.SETUP.Psi;
            u  = zeros(ssfem.SETUP.Ndof,ssfem.Rn);
            
            for j = 1:ssfem.SETUP.P
                psi_j = Psi{j}(ssfem.xi);
                d_j   = d(:,j);
                %                 u= u + d_j*psi_j';
                psi_j = repmat(psi_j',[ ssfem.SETUP.Ndof 1]);
                d_j = repmat(d_j,[1 ssfem.Rn]);
                u= u + d_j.*psi_j;
                printprogress(j,ssfem.SETUP.P);
            end
            
            ins=(1:ssfem.MeshData.NEdges);
            ins(ssfem.pecedge)=0;
            Et=zeros(ssfem.MeshData.NEdges,ssfem.Rn);
            Et(ins~=0,:)=u;
            ssfem.SETUP.Et = Et;
            
        end
        
     
        function trans = calcTrans_Stoch(obj,E)
            
            if(obj.evan)
                trans= 0;
                return;
            end
            
            OP_elems=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==obj.OP_PORT_NO).Port_elems];
            OP_face=[obj.MeshData.Faces(obj.MeshData.BoundaryFaces==obj.OP_PORT_NO).Port_faces];
            
            val=0;
            
            elems2nodes=obj.MeshData.Tetrahedron;
            elems2edges=obj.MeshData.Tet2Edges;
            nodes2coord=obj.MeshData.Nodes;
            faces2edges=obj.MeshData.Tri2Edges;
            faces2nodes=obj.MeshData.Triangle;
            signs=obj.MeshData.signs;
            T=obj.MeshData.Co;
            V=obj.MeshData.V;
            
            % function trans=scatp_mod(E,signs,nodes2coord,elems2nodes,elems2edges,belem_mark,rpedge,x_max,y_max,z_max,K)
            
            telems=OP_elems;
            
            
            for el=1:length(telems)
                ie=telems(el);
                fc=OP_face(el);
                
                % UL=U(:,:,:,:,ie);
                TL(:,:)=T(ie,:,:);
                
                
                % V(ie)=abs(V(ie));
                %%%%% write edge matrix
                %
                % edge=[A B; A C; A D; B C; D B; C D];
                local_node=[1 2; 1 3; 1 4; 2 3; 4 2; 3 4];
                
                
                fc_edge=faces2edges(fc,:);
                bedgeoftet=find(ismember(elems2edges(ie,:),fc_edge));
                fc_edge=elems2edges(ie,bedgeoftet);
                
                
                %     surf_nodes=unique(edges2nodes(elems2edges(ie,bedgeoftet),:));
                
                surf_nodes= faces2nodes(fc,:);
                SP=nodes2coord(surf_nodes,:);
                
                PQ=SP(2,:)-SP(1,:);
                PR=SP(3,:)-SP(1,:);
                n_hat=cross(PQ,PR)/norm(cross(PQ,PR));
                Pn=setdiff(elems2nodes(ie,:),surf_nodes);
                PS=nodes2coord(Pn,:)-SP(1,:);
                n_hat=n_hat*-sign(dot(n_hat,PS));
                
                
                sig=(signs(ie,bedgeoftet));
                
                
                
                i=bedgeoftet(1);j=bedgeoftet(2);k=bedgeoftet(3);
                % calculating local indices for nodes
                i1=local_node(i,1);
                i2=local_node(i,2);
                j1=local_node(j,1);
                j2=local_node(j,2);
                k1=local_node(k,1);
                k2=local_node(k,2);
                
                % calculating length of edges
                %                 li=1/pdist([edge(i,1:3);edge(i,4:6)],'euclidean');
                %                 lj=1/pdist([edge(j,1:3);edge(j,4:6)],'euclidean');
                %                 lk=1/pdist([edge(k,1:3);edge(k,4:6)],'euclidean');
                
                [li,lj,lk]=deal(1);
                
                %%%%%%%
                
                Li1=@(x,y,z)(1/(6*V(ie)))*(TL(i1,1)+x*TL(i1,2)+y*TL(i1,3)+z*TL(i1,4));
                Li2=@(x,y,z)(1/(6*V(ie)))*(TL(i2,1)+x*TL(i2,2)+y*TL(i2,3)+z*TL(i2,4));
                ci1=(1/(6*V(ie)))*TL(i1,2:4);%% [D E F]
                ci2=(1/(6*V(ie)))*TL(i2,2:4);%% [A B C]
                
                Lj1=@(x,y,z)(1/(6*V(ie)))*(TL(j1,1)+x*TL(j1,2)+y*TL(j1,3)+z*TL(j1,4));
                Lj2=@(x,y,z)(1/(6*V(ie)))*(TL(j2,1)+x*TL(j2,2)+y*TL(j2,3)+z*TL(j2,4));
                cj1=(1/(6*V(ie)))*TL(j1,2:4);%% [J K L]
                cj2=(1/(6*V(ie)))*TL(j2,2:4);%% [G H I]
                
                Lk1=@(x,y,z)(1/(6*V(ie)))*(TL(k1,1)+x*TL(k1,2)+y*TL(k1,3)+z*TL(k1,4));
                Lk2=@(x,y,z)(1/(6*V(ie)))*(TL(k2,1)+x*TL(k2,2)+y*TL(k2,3)+z*TL(k2,4));
                ck1=(1/(6*V(ie)))*TL(k1,2:4);%% [J K L]
                ck2=(1/(6*V(ie)))*TL(k2,2:4);%% [G H I]
                
                
                
                
                N1=@(x,y,z)li*(Li1(x,y,z)*ci2-Li2(x,y,z)*ci1);
                N2=@(x,y,z)lj*(Lj1(x,y,z)*cj2-Lj2(x,y,z)*cj1);
                N3=@(x,y,z)lk*(Lk1(x,y,z)*ck2-Lk2(x,y,z)*ck1);
                D1=@(x,y,z) -cross(n_hat,cross(n_hat,N1(x,y,z)));
                D2=@(x,y,z) -cross(n_hat,cross(n_hat,N2(x,y,z)));
                D3=@(x,y,z) -cross(n_hat,cross(n_hat,N3(x,y,z)));
                
                
                
                Ef =@(x,y,z) (sig(1)*repmat(E(fc_edge(1),:),3,1).*repmat(transpose(D1(x,y,z)),1,obj.Rn) +  ...
                    sig(2)*repmat(E(fc_edge(2),:),3,1).*repmat(transpose(D2(x,y,z)),1,obj.Rn) + ...
                    sig(3)*repmat(E(fc_edge(3),:),3,1).*repmat(transpose(D3(x,y,z)),1,obj.Rn));
                
                
                x=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),1);
                y=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),2);
                z=obj.MeshData.Nodes(obj.MeshData.Triangle(fc,:),3);
                z1=mean(z);
                Ef1=@(x,y) (sum(Ef(x,y,z1).^2,1));
                
                Ef1=@(x,y) vecnorm(Ef(x,y,z1),2).^2;
                
                in=obj.Ut.triquad(Ef1,4,x(1),x(2),x(3),y(1),y(2),y(3));
                
                val=val+in;
                
            end
            trans=val;
            
            trans=sqrt(trans/obj.ippow);
            %trans=abs(trans);
            
            
            
        end
        
        
        
    end
end