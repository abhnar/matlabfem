classdef KLE3D < handle
    properties
        MD
        Lambda
        Evec
        Bmat
        Cmat
        Cold
        Vol
        RAW
        Norms
        CORR
        KLDATA
        crfracs = [0.1 0.25 0.5 0.75 1 2 5 10];
        KLSet
        klcount = 0
    end
    methods
        function doms = getDomains(kle)

            for j=1:length(kle.KLSet)
                doms(j) = kle.KLSet{j}.domain;
            end
        end
        function sds = getSDs(kle)

            for j=1:length(kle.KLSet)
                sds(j) = kle.KLSet{j}.sd;
            end
        end
        function nkls = getNKLs(kle)

            for j=1:length(kle.KLSet)
                nkls(j) = kle.KLSet{j}.nkl;
            end
        end
        function N = getSdim(kle)
            N = 0;
            for j=1:length(kle.KLSet)
                N = N + kle.KLSet{j}.nkl;
            end
        end
        function ret = getKLEData(kle, corlenfrac)
            ret = [];
            ixs = zeros(length(kle.KLSet),1);
            if length(corlenfrac) == length(kle.KLSet)
                for i=1:length(corlenfrac)
                    ix  = find(kle.crfracs == corlenfrac);
                    if(isempty(ix))
                        error('Correlation length should obe 0.1, 0.25, 0.5, 0.75, 1,2,5 10.');
                    end
                    ixs(i) = ix;
                end
            elseif length(corlenfrac) == 1
                ix  = find(kle.crfracs == corlenfrac);
                if(isempty(ix))
                    error('Correlation length should obe 0.1, 0.25, 0.5, 0.75, 1,2,5 10.');
                end
                ixs = repmat(ix,size(kle.KLSet,2),1);
            end


            for i=1:length(kle.KLSet)
                kle.KLSet{i}.CORR = kle.KLSet{i}.CORR{ixs(i)};
                kle.KLSet{i}.corlenfactor = corlenfrac;
                ret =  kle.KLSet{i}.CORR;
            end

        end
        function setSD(kle, sds)
            if length(sds) ~= length(kle.KLSet)
                error(['There are ', num2str(length(kle.KLSet)), ' domains!']);
            end
            for i = 1:length(kle.KLSet)
                kle.KLSet{i}.sd = sds(i);
            end
        end
        function setNKL(kle, nkls)
            if length(nkls) ~= length(kle.KLSet)
                error(['There are ', num2str(length(kle.KLSet)), ' domains!']);
            end
            for i = 1:length(kle.KLSet)
                kle.KLSet{i}.nkl = nkls(i);
            end
        end
        function KL = process(kle, fem, domain_set)
            for domnum = 1:length(domain_set)
                domain = domain_set(domnum);
                tmp = split(fem.Meshfile,'.msh');
                klfname = strcat(tmp{1},'_KLE_',num2str(fem.MeshData.NT),'_',num2str(domain),'.mat');
                disp(['Attempting to access: "',klfname,'"'])
                if isfile(klfname)
                    L = load(klfname);
                    disp("Loaded from an existing KL evaluation");
                    KL = L.KL;
                    kle.KLSet{domnum} = KL;
                else
                    fprintf("Building KL Matrix system\n");
                    md = fem.MeshData;
                    ix1=find(md.TetType==domain);

                    tet=md.Tetrahedron(ix1,:);
                    V=md.V(ix1);
                    Co=md.Co(ix1,:,:);
                    t2e=md.Tet2Edges(ix1,:);
                    edges=t2e(:);
                    [edges,edix]=unique(edges(:));
                    [r,c]=ind2sub(size(t2e),edix);

                    ix = unique(tet(:));
                    val =1:length(ix);
                    map =sparse(ix,ix./ix,val);   %map_of_nodes

                    tet = map(tet);
                    nds = md.Nodes(ix,:);
                    xmx = max(nds(:,1));
                    ymx = max(nds(:,2));
                    zmx = max(nds(:,3));

                    xmn = min(nds(:,1));
                    ymn = min(nds(:,2));
                    zmn = min(nds(:,3));

                    %Construct passing structure
                    maxdim=max([xmx-xmn ymx-ymn zmx-zmn]);

                    passMD.Tetrahedron=tet;
                    passMD.Nodes=nds;
                    passMD.NT=size(tet,1);
                    passMD.NN=size(nds,1);
                    passMD.V=V;
                    passMD.Co=Co;
                    passMD.maxdim=maxdim;
                    kle.getKLE_multiple(passMD);
                    fprintf("Solving Eignen Value Problems\n");
                    kle.calcEig_Multiple(100);

                    supSet.V=V;
                    supSet.Co=Co;
                    supSet.NT=passMD.NT;
                    supSet.NTg=md.NT;
                    supSet.nds=nds;
                    supSet.NN = passMD.NN;
                    supSet.tet = tet;
                    supSet.edges = edges;
                    supSet.r = r;
                    supSet.ix1 = ix1;

                    DATA.supSet = supSet;
                    DATA.CORR = kle.CORR;

                    for it = 1:length(DATA.CORR)
                        %DATA.CORR{it}.domain = domain;
                        DATA.CORR{it}.corlen = maxdim*kle.crfracs(it);
                        DATA.CORR{it}
                    end

                    KL.supSet = DATA.supSet;
                    KL.CORR =  DATA.CORR;
                    KL.Volume = sum(V);
                    KL.domain = domain;
                    save(klfname,'KL');

                    kle.KLSet{domnum} = KL;
                end
                %              kle.evaluate_Phi();

            end
            KL = kle.KLSet;
        end
        function evaluate_Phi(kle, varargin)
            if(length(varargin) == 2)
                if(~strcmp(varargin{1},'nlambda'))
                    error("Input should be ''nlambda''");
                end
            end

            fprintf('Evaluating Eigen Function \n');
            for domnum = 1:length(kle.KLSet)
                a = kle.KLSet{domnum}.supSet.Co(:,:,1);
                b = kle.KLSet{domnum}.supSet.Co(:,:,2);
                c = kle.KLSet{domnum}.supSet.Co(:,:,3);
                d = kle.KLSet{domnum}.supSet.Co(:,:,4);


                ix1 = kle.KLSet{domnum}.supSet.ix1;
                nds = kle.KLSet{domnum}.supSet.nds;
                tet = kle.KLSet{domnum}.supSet.tet;
                V = kle.KLSet{domnum}.supSet.V;

                for it = 1:length(kle.KLSet{domnum}.CORR)
                    %                     normset = kle.KLSet{domnum}.CORR{it}.Norms;
                    llim = min(length(kle.KLSet{domnum}.CORR{it}.Lambda), varargin{2});
                    fprintf('for fraction %f\n',kle.crfracs(it));
                    Phi = {};
                    for j = 1: llim
                        eigvec = kle.KLSet{domnum}.CORR{it}.Evec(:,j);

                        phi_centroid=zeros(kle.KLSet{domnum}.supSet.NTg,1);
                        val = 0;
                        for i=1:length(tet)
                            gix = ix1(i);
                            lix = i;
                            Xc=mean(nds(tet(lix,:),:));

                            field_centroid=0;
                            for k=1:4
                                Lk = @(x) 1/(6*V(lix))*(a(lix,k)+b(lix,k)*x(1) +c(lix,k)*x(2)+d(lix,k)*x(3));

                                field_centroid = field_centroid+Lk(Xc)*eigvec(tet(lix,k));
                            end
                            integ=0;
                            for k=1:4
                                Lk = @(x) 1/(6*V(lix))*(a(lix,k)+b(lix,k)*x(1) +c(lix,k)*x(2)+d(lix,k)*x(3));

                                integ = integ+Lk(Xc)*eigvec(tet(lix,k))*eigvec(tet(lix,k));
                            end


                            val=val+integ*V(lix);
                            %disp([edges(i) md.Tet2Edges(ix1(r(i)),:)])

                            phi_centroid(gix)=field_centroid;

                            printprogress(i + (j-1)*size(tet,1),size(tet,1)*llim);
                        end
                        Phi{j} = phi_centroid/sqrt(val);
                        %                         disp(sqrt(val));
                    end
                    kle.KLSet{domnum}.CORR{it}.Phi =  Phi;
                end
            end

        end
        function getKLE_multiple(obj,mshdata)
            %             ME=MeshEngine2;
            %             obj.MD=ME.readNodalMesh(fname);
            obj.MD=mshdata;
            md=obj.MD;
            nodes=md.Nodes;
            tet=md.Tetrahedron;
            Ut=Utilities;
            V=md.V;
            a=md.Co(:,:,1);
            b=md.Co(:,:,2);
            c=md.Co(:,:,3);
            d=md.Co(:,:,4);
            obj.Vol=sum(V);

            B = zeros(md.NN);
            [points,W]=Ut.tetweights(4,'keast');
            NP=length(W);
            NT=length(tet);


            xmx=max(nodes(:,1));
            ymx=max(nodes(:,2));
            zmx=max(nodes(:,3));

            xmn=min(nodes(:,1));
            ymn=min(nodes(:,2));
            zmn=min(nodes(:,3));

            maxdim=max([xmx-xmn ymx-ymn zmx-zmn]);


            for e=1:md.NT
                Be=zeros(4,4);
                vertex=[nodes(tet(e,:),1)  nodes(tet(e,:),2)  nodes(tet(e,:),3)];
                %vertex=[nodes(tet(e,1),:)  nodes(tet(e,2),:)  nodes(tet(e,3),:)];
                for i=1:4
                    Li=@(x) 1/(6*V(e))*(a(e,i)+b(e,i)*x(1) +c(e,i)*x(2)+d(e,i)*x(3));


                    for j=1:4
                        Lj=@(x) 1/(6*V(e))*(a(e,j)+b(e,j)*x(1) +c(e,j)*x(2)+d(e,j)*x(3));


                        fun = @(x) Li(x)*Lj(x);
                        val=0;
                        for p = 1:NP
                            X=points(p,:);

                            X1=X(1)*vertex(1,:)+X(2)*vertex(2,:)+X(3)*vertex(3,:)+X(4)*vertex(4,:);
                            val = val + fun(X1)*W(p);
                        end
                        Be(i,j)=val*V(e);

                    end
                end

                B(tet(e,:),tet(e,:))=B(tet(e,:),tet(e,:))+Be;
            end









            covlength=maxdim;
            fprintf("Covariance length: %d\n",covlength);
            cov =@(x1,x2) exp(-vecnorm(x1-x2)/(covlength));
            %cov =@(x1,x2) ((x1(2)<0.5 && x2(2)<0.5) || (x1(2)>=0.5 && x2(2)>=0.5));


            vale = zeros(4,NT,NP);
            coordmat = zeros(NT,3,NP);
            for e = 1:NT
                vertex=[nodes(tet(e,:),1)  nodes(tet(e,:),2)  nodes(tet(e,:),3)];
                for j=1:4
                    Lj=@(x) 1/(6*V(e))*(a(e,j)+b(e,j)*x(1) +c(e,j)*x(2)+d(e,j)*x(3));
                    for p =1:NP
                        X=points(p,:);
                        X1=X(1)*vertex(1,:)+X(2)*vertex(2,:)+X(3)*vertex(3,:)+X(4)*vertex(4,:);
                        vale(j,e,p)=Lj(X1)*W(p)*V(e);
                    end
                end

                for p =1:NP
                    X=points(p,:);
                    X1=X(1)*vertex(1,:)+X(2)*vertex(2,:)+X(3)*vertex(3,:)+X(4)*vertex(4,:);
                    coordmat(e,:,p)=X1;
                end

            end










            obj.Bmat = B;
            obj.Cmat = buildKLMatrices(vale,coordmat,NT,NP,md.NN,tet,covlength);



        end
        function ret =calcEig_Multiple(obj,n)
            obj.CORR = {};
            B=obj.Bmat;
            for it = 1:length(obj.Cmat)
                fprintf('For corrlen frac %2.2f\n',obj.crfracs(it));
                C=obj.Cmat{it};
                if(size(B,1)<n)
                    fprintf('N should be less than %d. Setting n= %d.\n',size(B,1),size(B,1));
                    n=size(B,1);
                end
                [eigvec,eigval]=eigs(C,B,size(B,1));
                obj.RAW.eigvec=eigvec;
                obj.RAW.eigval=eigval;

                [eval1,ix]=sort(diag(eigval),'desc');


                obj.Lambda = eval1;
                obj.Evec = eigvec(:,ix);

                obj.CORR{it}.Lambda=eval1;
                obj.CORR{it}.Evec=eigvec(:,ix);

                ret=eval1(1:n);


            end
        end
        function val= checkorthogonality(obj,i1,i2)
            MD = obj.MD;
            tets = MD.Tetrahedron;
            nds = MD.Nodes;
            NT = MD.NT;

            Co = MD.Co;
            a = Co(:,:,1);
            b = Co(:,:,2);
            c = Co(:,:,3);
            d = Co(:,:,4);

            V = MD.V;
            evec = obj.Evec;
            eval = obj.Lambda;

            val=0;
            for e=1:NT
                field=0;
                field1 = 0;
                field2 = 0;
                for k=1:4
                    Lk = @(x) 1/(6*V(e))*(a(e,k)+b(e,k)*x(1) +c(e,k)*x(2)+d(e,k)*x(3));
                    xyz=mean(nds(tets(e,:),:));
                    field = field+Lk(xyz)*evec(tets(e,k),i1)*evec(tets(e,k),i2)*Lk(xyz);

                end

                val=val+field*V(e);


            end

            %h=figure;
            %disp(val);

        end
    end
end

