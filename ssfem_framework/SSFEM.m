classdef SSFEM<  SFEM
    properties

        d_stoch;
        show_progress = true;

        
    end
   
    methods
        
        
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
        
        function ssfemkle(ssfem,f)
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
            clear X Y VALS S x y val;
            
            G=zeros(Ndof*P,1);
            G(1:Ndof)=b;
            fprintf('Solving\n');
            tic;
            xr = ST\G;
            t1 = toc;
            cprintf('*black','\tSolution time     : %f sec\n', t1);
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
            ssfem.RESULTS.N = size(ST,1);
            ssfem.RESULTS.NZ = nnz(ST);
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
        
        function ssfemrv(ssfem,f)
            
            ssfem.buildSystem();
            [K, b]=ssfem.buildDeterministicSystem(f);
            T =-ssfem.K0^2* ssfem.RV.T;
          
            Tset = cell(2,1);
            Tset{1}=K;
            Tset{2}=T;
            [cijk,P]=PC.c_ijk(1,ssfem.p_order,1);
            ssfem.RESULTS.P = ssfem.p_order;
            ssfem.RESULTS.type = 'RV';
            Ndof = size(K,1);
            Stemp = cell(P,P);
            ST = cell(P);
            for j = 1 : P
                
                for k = 1 : P
                    S = sparse(Ndof,Ndof);
                    
                    for i =1 : 2
                        S =  S + cijk{i}(j,k) * Tset{i};
                    end
                    Stemp{j,k}=S;
                    
                    
                end
                
                ST{j}=horzcat(Stemp{j,:});
                
            end

            ST=vertcat(ST{:});
            clear Stemp;
            
            G=zeros(Ndof*P,1);
            G(1:Ndof)=b;
            disp('Solving');
            tic;
            xr = ST\G;
            t1 = toc;
            ssfem.d_stoch = reshape(xr,Ndof,P);
            ssfem.RESULTS.solvetime = t1;
            ssfem.RESULTS.N = size(ST,1);
            ssfem.RESULTS.NZ = nnz(ST);
            d = ssfem.d_stoch;
            disp('Evaluating PCE');
            Psi=PC.getHermite_PC(1,ssfem.p_order);
            u  = zeros(Ndof,ssfem.Rn);
            
            for j = 1:P
                psi_j = Psi{j}(ssfem.xi);
                d_j   = d(:,j);
                u= u + d_j*psi_j';
                printprogress(j,P);
            end
            
            ins=(1:ssfem.MeshData.NEdges);
            ins(ssfem.pecedge)=0;
            Et=zeros(ssfem.MeshData.NEdges,ssfem.Rn);
            Et(ins~=0,:)=u;
            
            ssfem.calcTrans_Stoch(Et);
            ssfem.RESULTS.f = f;
            ssfem.RESULTS.seed = ssfem.seed;
            
            ssfem.pushResult();
            
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