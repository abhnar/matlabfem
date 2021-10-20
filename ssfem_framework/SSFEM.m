classdef SSFEM<  SFEM
    properties

        d_stoch;

        
    end
   
    methods
        
        
        function assignRandomMaterialVariation(ssfem, domain, sd)
            
            assignRandomMaterialVariation@SFEM(ssfem, domain, sd, 'SSFEM');
            
            ssfem.buildSystem();
            
            epr=ssfem.MeshData.TetType*0;
            
            temp=ssfem.MeshData.TetType;
            temp(temp~=domain)=0;
            temp(temp==domain)=1;
            epr = sd.*(temp);
            
            ssfem.Assemble(epr);
            
            ssfem.RV.T = ssfem.T;
            ssfem.RV.M = ssfem.M;
            ssfem.RV.T(ssfem.pecedge,:) = [];
            ssfem.RV.T(:,ssfem.pecedge) = [];
            ssfem.RV.M(ssfem.pecedge,:) = [];
            ssfem.RV.M(:,ssfem.pecedge) = [];
            
        end
        
        function assignSpatialMaterialVariation(ssfem, kle)
            
              
            assignSpatialMaterialVariation@SFEM(ssfem, kle, 'SSFEM');
            
            
            
            ssfem.buildSystem();
            
            
            
            for i=2:ssfem.nkl+1
                norm_i = ssfem.getNorm(kle,i-1);
                epr = kle.KLDATA.sd*sqrt(kle.CORR.Lambda(i-1))*kle.CORR.Phi{i-1}/norm_i;
                ssfem.Assemble(epr);
                ssfem.KLE.Tset{i} =  ssfem.T;
                ssfem.KLE.Tset{i}(:,ssfem.pecedge)=[];
                ssfem.KLE.Tset{i}(ssfem.pecedge,:)=[];
            end
            
            
            
            
            
        end
        
        function ssfemkle(ssfem,f)
            ssfem.buildSystem();
            [K, b]=ssfem.buildDeterministicSystem(f);
            
            Tset = ssfem.KLE.Tset;
            Tset{1} = K;
            for i=2:length(Tset)
                Tset{i} = -ssfem.K0^2*Tset{i};
            end
            [cijk,P]=PC.c_ijk(ssfem.nkl,ssfem.p_order,1);
            ssfem.RESULTS.P=P;
            ssfem.RESULTS.p_order=ssfem.p_order;
            ssfem.RESULTS.KLterms=ssfem.nkl;
            ssfem.RESULTS.type = 'KLE';
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
            
            ssfem.d_stoch = reshape(xr,Ndof,P);
            ssfem.RESULTS.solvetime = t1;
            ssfem.RESULTS.N = size(ST,1);
            ssfem.RESULTS.NZ = nnz(ST);
            d = ssfem.d_stoch;
            disp('Evaluating PCE');
            Psi=PC.getHermite_PC(ssfem.nkl,ssfem.p_order);
            u  = zeros(Ndof,ssfem.Rn);
            rng(ssfem.seed);
            ssfem.xi = randn(ssfem.Rn,ssfem.nkl);
            for j = 1:P
                psi_j = Psi{j}(ssfem.xi);
                d_j   = d(:,j);
                %                 u= u + d_j*psi_j';
                psi_j = repmat(psi_j',[ Ndof 1]);
                d_j = repmat(d_j,[1 ssfem.Rn]);
                u= u + d_j.*psi_j;
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
        
        function ssfemrv(ssfem,f)
            
            ssfem.buildSystem();
            [K, b]=ssfem.buildDeterministicSystem(f);
            T =-ssfem.K0^2* ssfem.RV.T;
            M = ssfem.RV.M;
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
            rng(ssfem.seed);
            ssfem.xi = randn(ssfem.Rn,1);
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
            obj.RESULTS.T=trans;
            
            
        end
        
        
        
    end
end