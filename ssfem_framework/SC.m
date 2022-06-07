classdef SC<handle
    %Author: Abhijith B N
    %email: abhijithbn@mail.com
    properties
        SDIM = 1;
        MN= 0;
        SD = 1;
        GRID;
        NSC = 0;
        Rn = 20000;
        ncpt = 0;
        xisc = 0;
        indices = [];
        GRID_NORMALIZED = [];

    end
    methods
        function setSDIM(obj, sdim)
            obj.SDIM = sdim;
        end
        function setMeans(obj, mns)
            obj.MN = mns;
        end
        function setSDs(obj, sds)
            obj.SD = sds;
        end
        function grid = generateCollocationPoints(obj, ncpt)
            if(length(ncpt)==1)
                ncpt = repmat(ncpt,length(obj.MN),1);
            end
            if(length(obj.SD)~=length(obj.MN))
                error('Different number of means and stds');
            end
            obj.SDIM=length(obj.MN);
            xi = cell( obj.SDIM,1);
            obj.xisc = xi;

            for i=1:obj.SDIM
                obj.xisc{i} = obj.hermitequads(ncpt(i));
                xi{i} = obj.MN(i)+ obj.xisc{i}*obj.SD(i);

            end
            grid = obj.cartprod(xi);
            obj.GRID_NORMALIZED = obj.cartprod( obj.xisc);

            obj.GRID = grid;
            obj.NSC  = size(grid,1);
            obj.indices = zeros(size(grid));
            for k=1:obj.NSC
                for i=1:obj.SDIM
                    obj.indices(k,i) =  find(xi{i} == grid(k,i)) ;
                end
            end
            obj.ncpt = ncpt;

        end

        function y0 = collocate(obj, y_sc_sampled)
            xi = cell(obj.SDIM,1);
            rng(1);
            for i=1:obj.SDIM
                xi{i} = randn(obj.Rn,1);
            end

            y0 = zeros(length(xi{1}),1);

            for k=1:obj.NSC

                pr = 1;
                for i=1:obj.SDIM
                    j=obj.indices(k,i);
                    pr = pr.*obj.lagrange_basis_function_1d(obj.ncpt(i)-1,obj.xisc{i},j,xi{i});
                end
                y0 = y0 + y_sc_sampled(k)*pr;

            end


        end
        function xp=hermitequads(obj,K)
            [~,H,~,~,~]=obj.Hermite_PC(1,K+1);
            A=((solve(H{K+1})));
            xp=eval(vpa(A));
            xp = sort(xp);
        end
        function yi = lagrange_basis_function_1d (obj, mx, xd, i, xi )

            %*****************************************************************************80
            %
            %% LAGRANGE_BASIS_FUNCTION_1D evaluates a 1D Lagrange basis function.
            %
            %  Licensing:
            %
            %    This code is distributed under the GNU LGPL license.
            %
            %  Modified:
            %
            %    04 August 2012
            %
            %  Author:
            %
            %    John Burkardt
            %
            %  Input:
            %
            %    integer MX, the degree of the basis function.
            %
            %    real XD(MX+1), the interpolation nodes.
            %
            %    integer I, the index of the basis function.
            %    1 <= I <= MX+1.
            %
            %    real XI, the evaluation point.
            %
            %  Output:
            %
            %    real YI, the value of the I-th Lagrange 1D basis function
            %    for the nodes XD, evaluated at XI.
            %
            if(size(xd,1)>size(xd,2))
                xd = xd';
            end
            if ( xi == xd(i) )
                yi = 1.0;
            else
                j = find ( 1 : mx + 1 ~= i );
                %yi = prod ( xi - xd(j) ) / prod ( xd(i) - xd(j) );
                yi = prod(repmat(xi,[1 size(j,2)]) - xd(j),2)/prod ( xd(i) - xd(j) );
            end

            return
        end

        function X = cartprod(obj, varargin)

            varargin = varargin{1};
            numSets = length(varargin);
            for i = 1:numSets,
                thisSet = sort(varargin{i});
                if ~isequal(prod(size(thisSet)),length(thisSet)),
                    error('All inputs must be vectors.')
                end
                if ~isnumeric(thisSet),
                    error('All inputs must be numeric.')
                end
                if ~isequal(thisSet,unique(thisSet)),
                    error(['Input set' ' ' num2str(i) ' ' 'contains duplicated elements.'])
                end
                sizeThisSet(i) = length(thisSet);
                varargin{i} = thisSet;
            end
            X = zeros(prod(sizeThisSet),numSets);
            for i = 1:size(X,1),

                % Envision imaginary n-d array with dimension "sizeThisSet" ...
                % = length(varargin{1}) x length(varargin{2}) x ...

                ixVect = obj.ind2subVect(sizeThisSet,i);

                for j = 1:numSets,
                    X(i,j) = varargin{j}(ixVect(j));
                end
            end



        end

        function X = ind2subVect(obj, siz,ndx)

            n = length(siz);
            k = [1 cumprod(siz(1:end-1))];
            ndx = ndx - 1;
            for i = n:-1:1,
                X(i) = floor(ndx/k(i))+1;      % replaced "varargout{i}" with "X(i)"
                ndx = rem(ndx,k(i));
            end
        end

        function [alpha,Psi_s,Psi_p,PsiSqNorm,P] = Hermite_PC(obj,M,p_order)

            %% Calculate the basis size of Psi
            P = 1;
            for s = 1:p_order
                P = P + (1/factorial(s))*prod(M+(0:s-1));   % Eq.(3.13) Ref.(2)
            end

            %% Calculate 1D Hermite polynomials: Recurrence relation
            % symbolic
            syms xi;
            He_s    = cell(p_order,1);
            He_s{1} = sym(1);
            He_s{2} = xi;
            for j = 2:p_order+1
                He_s{j+1} = expand(xi*He_s{j} - (j-1)*He_s{j-1});
            end
            % polynomial
            He_p    = cell(p_order,1);
            He_p{1} = 1;       % H_1 = 1
            He_p{2} = [1 0];   % H_2 = x
            for n = 2:p_order+1
                He_p{n+1} = [He_p{n} 0] - (n-1)*[0 0 He_p{n-1}];   % recursive formula
            end

            %% Define the number of RVs
            x   = cell(1,M);
            H_s = cell(p_order,M);   % Hermite polinomial for each dimension
            H_p = cell(p_order,M);   % Hermite polinomial for each dimension
            for j = 1:M
                x{j} = sym(sprintf('xi_%d',j));
                for i = 1:p_order+1
                    H_s{i,j} = subs(He_s{i},xi,x{j});
                    H_p{i,j} = He_p{i};
                end
            end

            %% M-dimensional PC computation
            Psi_s  = cell(P,1);   % symbolic version
            Psi_p  = cell(P,1);   % polynomial version
            alpha  = PC.multi_index(M,p_order);  % create the multi-index
            for i = 2:P+1
                mult_s = 1;
                mult_p = 1;
                for j = 1:M
                    mult_s = mult_s*H_s{alpha(i-1,j)+1,j};
                    mult_p = conv(mult_p,H_p{alpha(i-1,j)+1,j});
                end
                Psi_s{i-1} = mult_s;
                Psi_p{i-1} = mult_p;
            end

            %% Calculate the square norm
            PsiSqNorm  = prod(factorial(alpha),2);

            %% show results
            %fprintf('Number of random variables (K-L terms): M = %d \n',M);
            %fprintf('Order of the polynomial chaos: p = %d \n',p_order);
            %fprintf('Total number of polynomial chaos: p = %d \n\n',P);
            for k = 1:P
                %fprintf('j = %d \t',k-1);
                % fprintf('E[Psi^2] = %d \t\t',PsiSqNorm(k));
                %fprintf('Psi   = %s \n',char(Psi_s{k}));
            end

            return;
        end
    end
end

function alpha = multi_index(M,p)
%% multi-index sequence for the computation of M-dimensional polynomials
%{
    --------------------------------------------------------------------------
    Created by:                       Date:           Comment:
    Felipe Uribe                      Oct/2014        ---
    furibec@unal.edu.co
    Universidad Nacional de Colombia
    Manizales Campus
    --------------------------------------------------------------------------
    *Input:
     M           % number of K-L terms (number of random variables)
     p_order     % order of PC
    --------------------------------------------------------------------------
    *Output:
     alpha       % multi-index sequence
    --------------------------------------------------------------------------
    Based on:
    1."Numerical methods for stochastic computations. A spectral method approach"
       D. Xiu. 2010. Princeton University Press.
    2."Stochastic finite element methods and reliability"
       B. Sudret and A. Der Kiureghian. State of the art report.
    --------------------------------------------------------------------------
%}

%% procedure
alpha    = cell(p+1,1);    % multi-index
alpha{1} = zeros(1,M);     % multi-index for length 0

switch M;
    case 1;  % dimension = 1
        for q = 1:p
            alpha{q+1} = q;
        end
    otherwise;  % dimension>1
        for q = 1:p
            s          = nchoosek(1:M+q-1,M-1);
            s1         = zeros(size(s,1),1);
            s2         = (M+q)+s1;
            %disp(s)
            alpha{q+1} = flipud(diff([s1 s s2],1,2))-1;   % -1 due to MATLAB indexing
            if sum(alpha{q+1},2) ~= q*ones(nchoosek(M+q-1,M-1),1)
                error('The sum of each row has to be equal to q-th order');
            end
        end
end
alpha = cell2mat(alpha);   % as in Ref.(1)
% alpha = flipud(alpha')';   % ~as in Ref.(2)

return;
end


