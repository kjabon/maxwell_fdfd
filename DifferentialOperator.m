classdef DifferentialOperator < handle
    %DIFFERENTIALOPERATOR Summary of this class goes here
    %   Later: Also pass in desired dimension (1: d/dx, 2: d/dy, etc.)
    %   Also pass in originSubGrid, i.e. whether we start on 's' or 't'
    %   (star/sircle, or triangle) 
    
    properties
        domain %dx, NxS, NxT, x, BC, numDimensions
        dxOriginSubgrid
        dyOriginSubgrid
        dzOriginSubgrid
        dxst %d/dx operator going from s to t grid
        dxts %d/dx operator going from t to s grid
        dyst %d/dy operator going from s to t grid
        dyts %d/dy operator going from t to s grid
        dzst %d/dz operator going from s to t grid
        dzts %d/dz operator going from t to s grid
        Sxst %Sx (PML) matrix going from s to t grid 
        Sxts %Sx (PML) matrix going from t to s grid
        Syst %Sy (PML) matrix going from s to t grid
        Syts %Sy (PML) matrix going from t to s grid
        Szst %Sz (PML) matrix going from s to t grid
        Szts %Sz (PML) matrix going from t to s grid
        
    end
    
    methods
        function obj = DifferentialOperator(domain)
            %UNTITLED3 Construct an instance of this class
            %   Initial call returns first derivative, with specified
            %   originSubgrid
            obj.domain = domain;
            
            %create and save all operators
            obj.dxOriginSubgrid = 'uninitialized';
            obj.dyOriginSubgrid = 'uninitialized';
            obj.dzOriginSubgrid = 'uninitialized';
            
            %Assume we have at least one dimension...
            obj.dxts = obj.dx('t');
            obj.dxst = obj.dx('s');
            
            if domain.numDimensions >= 2
                %initialize dy
                obj.dyts = obj.dy('t');
                obj.dyst = obj.dy('s');
             
            end
            
            %later: optionally only create up to n dimensions, or just
            %create on demand when called
        end
        function laplacian = Laplacian(obj)
            dxSquared = obj.dxts*obj.dxst;
            dySquared = obj.dyts*obj.dyst;
            laplacian = dxSquared+dySquared; 
        end
        function operatorMatrix = dx(obj,originSubgrid)
            
            %METHOD1 Returns d/dx of current operator matrix
            %   Detailed explanation goes here
            domain = obj.domain;
            
            if domain.numDimensions == 1
                operatorMatrix = obj.dx1D(originSubgrid);
            elseif domain.numDimensions == 2
                operatorMatrix = obj.dx2D(originSubgrid);
            end
            
        end
        
        %These should be Private if possible in matlab
        function operatorMatrix = dx1D(obj, originSubgrid)
             %S (PML) matrix going from s to t grid
            domain = obj.domain;  %#ok<*PROPLC,*PROP>
            BC = domain.BC;
            x = domain.x;
            
            %% Handle arguments
            if nargin == 1 %originSubgrid arg is blank
                if ~(strcmp(obj.dxOriginSubgrid,'s') || strcmp(obj.dxOriginSubgrid,'t'))
                    error('OriginSubgrid has not yet been set to s or t.');
                end

                %swap the originSubgrid
                if obj.dxOriginSubgrid == 's'
                    obj.dxOriginSubgrid = 't'; 
                elseif obj.dxOriginSubgrid == 't'
                    obj.dxOriginSubgrid = 's'; 
                end
            elseif originSubgrid == 0
                if ~(strcmp(obj.dxOriginSubgrid,'s') || strcmp(obj.dxOriginSubgrid,'t'))
                    error('OriginSubgrid has not yet been set to s or t.');
                end
                %do not swap
            elseif strcmpi(originSubgrid, 's') || strcmpi(originSubgrid, 'e')
                obj.dxOriginSubgrid = 's'; 
                originSubgrid = 's';
                if ~isempty(obj.dxst)
                    operatorMatrix = obj.dxst;
                    return;
                end
                %otherwise, make the matrix and return it
            elseif strcmpi(originSubgrid, 't') || strcmpi(originSubgrid, 'h')
                obj.dxOriginSubgrid = 't'; 
                originSubgrid = 't';
                if ~isempty(obj.dxts)
                    operatorMatrix = obj.dxts;
                    return;
                end
            end
            
            %% Construct components for sparse matrix

            %i is first dimension, or row # (also "destination dimension" 
            %of differential matrix)
            %j is second dimension, or column # (also "origin dimension" of
            %differential matrix)
            if strcmpi(obj.dxOriginSubgrid, 's')
                s = true;
            else 
                s = false;
            end
                t = ~s;

            if s
                Ni = domain.NxT(1);
                Nj = domain.NxS(1);
            else %t
                Ni = domain.NxS(1);
                Nj = domain.NxT(1);
            end

            %Set up indices for sparse matrix
            %(nothing intuitive, extracted from writing matrices by hand)
            if ((BC(1,1) == BCs.antiSymS || BC(1,1) == BCs.symT) && s )||...    %PEC
                ((BC(1,1) == BCs.antiSymT || BC(1,1) == BCs.symS) && t )        %PMC
                i11 = 1;
                j11 = 1;
                i21 = 2;
                j21 = 1;
            elseif ((BC(1,1) == BCs.symS || BC(1,1) == BCs.antiSymT) && s )||...%PMC
                    ((BC(1,1) == BCs.symT || BC(1,1) == BCs.antiSymS) && t )    %PEC
                i11 = 1;
                j11 = 2;
                i21 = 1;
                j21 = 1;
            end

            if ((BC(1,2) == BCs.antiSymS || BC(1,2) == BCs.symT) && s)||...     %PEC 
                    ((BC(1,2) == BCs.antiSymT || BC(1,2) == BCs.symS) && t)     %PMC
                i1end = Ni-1;
                j1end = Nj;
                i2end = Ni;
                j2end = Nj;
            elseif ((BC(1,2) == BCs.symS || BC(1,2) == BCs.antiSymT)&&s)||...   %PMC
                    ((BC(1,2) == BCs.symT || BC(1,2) == BCs.antiSymS)&&t)       %PEC
                i1end = Ni;
                j1end = Nj;
                i2end = Ni;
                j2end = Nj-1; 
            end


            i1 = i11:i1end;
            j1 = j11:j1end;
            i2 = i21:i2end;
            j2 = j21:j2end;

            %Build the first diagonal for this matrix
            vec1 = ones(size(i1)); 
            %Build the second diagonal for this matrix
            vec2 = -1*ones(size(i2));  

            %% Cases for left boundary (editing diagonals):

%             if BC(1) == BCs.antiSymS 
%                 vec2(1) = 0;        
%             elseif BC(1) == BCs.antiSymT 
%                 vec1(1) = 0;
            if BC(1,1) == BCs.symS 
                if s
                    vec2(1) = sqrt(2)*vec2(1);
                else
                    vec1(1) = sqrt(2)*vec1(1);
                end
                %vec2(1) = -sqrt(2);
            elseif  BC(1,1) == BCs.symT 
                if s
                    vec1(1) = vec1(1)*sqrt(2);
                else
                    vec2(1) = vec2(1)*sqrt(2);
                end
                %vec1(1) = sqrt(2);    
            end

            %% Cases for right boundary (editing diagonals):

%             if BC(2) == BCs.antiSymS 
%                 vec1(end) = 0;
%             elseif  BC(2) == BCs.antiSymT 
%                 vec2(end) = 0;
%             else
            if BC(1,2) == BCs.symS
                if s
                    vec1(end) = sqrt(2)*vec1(end);
                else
                    vec2(end) = sqrt(2)*vec2(end);
                end
                %vec1(end) = sqrt(2);
            elseif  BC(1,2) == BCs.symT 
                if s
                    vec2(end) = sqrt(2)*vec2(end);
                else
                    vec1(end) = sqrt(2)*vec1(end);
                end
                %vec2(end) = -sqrt(2);
            end



            %% Combine left and right BCs into one matrix


            i = [i1,i2];
            j = [j1,j2];
            combinedVec = [vec1,vec2]'; %Create one long vector

            %% But first! Modify matrix based on local dx

             %Index of first s point in x vector
            if BC(1,1) == BCs.antiSymT || BC(1,1) == BCs.symS %Start on s point in x direction
                s1 = 1; %Not used, we use -op'
                t1 = 2;
            else %Start on t point in x direction
                s1 = 2; %Not used, we use -op'
                t1 = 1;
            end

            %Build dx vector
            %This is a little difficult to read, but was necessary in
            %order to use no for-loops

            if s
                x1 = x(s1:2:(Nj-1)*2); 
                x2 = x(s1+2:2:(Nj)*2); %Note the end index here may overshoot by 1, but matlab won't go the extra half cell because of the :2:


                if BC(1,1) == BCs.antiSymT || BC(1,1) == BCs.symS %End on an s point
                    dx1 = [];
                else %End on a t point
                    dx1 = (x(2) - x(1))*2;
                end


                if BC(1,2) == BCs.antiSymT || BC(1,2) == BCs.symS %End on an s point
                    dxend = [];
                else %End on a t point
                    dxend = (x(end) - x(end-1))*2; 
                end

            else
                x1 = x(t1:2:(Nj-1)*2); 
                x2 = x(t1+2:2:(Nj)*2); %Note the end index here may overshoot by 1, but matlab won't go the extra half cell because of the :2:

                if BC(1,1) == BCs.antiSymT || BC(1,1) == BCs.symS %End on an s point
                    dx1 = (x(2) - x(1))*2;
                else %End on a t point
                    dx1 = [];
                end

                if BC(1,2) == BCs.antiSymT || BC(1,2) == BCs.symS %End on an s point
                    dxend = (x(end) - x(end-1))*2; 
                else %End on an t point
                    dxend = [];
                end
            end

            dx = x2-x1; %contains the dx from surrounding originSubgrid points, around target points
            dx = [dx1 dx dxend];
            
%% Old PML addition (with for loop)
% %%%%% Old, more legible way of computing dx with for loops %%%%%%%%%%%%%
%                 %We can compute a dx vector, which corresponds to the dx = x(i+1)-x(i-1)
%                 for ii = 2:length(x)-1
%                    dxold(ii-1) = x(ii+1)-x(ii-1); 
%                 end
%                 
%                 dx1 = (x(2) - x(1))*2;
%                 dxend = (x(end) - x(end-1))*2;
%                 dxold = [dx1 dxold dxend];
%                 %Our dx at the two endpoints is always (x(2)-x(1))*2
% 
%                 %We must access this correctly, because this corresponds to
%                 %both s and t points
% 
%                 %So we can calculate the whole thing, then create dxst (which
%                 %takes the dx between neighboring s points), with dxst =
%                 %dx(t1:2:NxT), or %dxts = %dx(s1:2:NxS)
%                 if s
%                     dxold = dxold(t1:2:Ni*2); %origin s, destination t
%                 else
%                     dxold = dxold(s1:2:Ni*2); %origin t, destination s
%                 end
% 
% %%%%%%%%%%%      %Then we do operatorMatrix = sparse(...)./dxst' %%%%%%%%

%% Add PMLs in x vector via complex coordinate stretching
            power = 3;
            %^Later, make this a param

            Sx = ones(size(domain.x));

            if domain.PML(1,1) > 0 %Left boundary 
                % TODO: Should increase as s(x)/omega so all frequencies will
                % decay at same rate
                % s(x) should be a cubic function for starters
                % Increase up to max value: PMLMag

                imComponent = (domain.PML(1,1)*domain.dx(1)):-domain.dx(1)/2:dx(1)/2;    % E.g.: 20:0.5:0.5
                A = domain.PMLMag /((domain.PML(1,1)*domain.dx(1))^power);               % Calculate scaling factor to end at PMLMag 
                imComponent = 1i*(A*imComponent.^power);                    % Generate cubically-increasing imaginary component
                %imComponent = -1* imComponent;
                %Instead of adding imaginary component directly to dx
                %vector, encapsulate it in a Sst diagonal matrix (and
                %Sts), then later selectively front multiply (Dxst = Sst*Dxst)

                Sx(1:domain.PML(1,1)*2) = Sx(1:domain.PML(1,1)*2) + imComponent;
                %domain.x(1:domain.PML(1,1)*2) = domain.x(1:domain.PML(1,1)*2) + imComponent; 
                                                                        % ^Overlay generated complex coordinates onto
                                                                        % left side of existing x vector
            end

            if domain.PML(1,2) > 0 %Right boundary
                imComponent = domain.dx(1)/2:domain.dx(1)/2:domain.PML(1,2)*domain.dx(1);
                A = domain.PMLMag /((domain.PML(1,2)*domain.dx(1))^power);
                imComponent = 1i*(A*imComponent.^power);
                %imComponent = -1* imComponent;

                Sx(end-domain.PML(1,2)*2+1:end) = Sx(end-domain.PML(1,2)*2+1:end) + imComponent;

                %domain.x(end-domain.PML(1,2)*2+1:end) = domain.x(end-domain.PML(1,2)*2+1:end) + imComponent;
            end
            
            %Get Sx for DESTINATION subgrid (not entirely certain why, but
            %the dimensions work)
            if s
                Sx = domain.GetSubvector('t', 1, Sx);
            else
                Sx = domain.GetSubvector('s', 1, Sx);
            end
            
            Sx = diag(sparse((1./Sx)));
            %Do not need to replicate into 2d; this is taken care of by calling function
       
            if s
                obj.Sxst = Sx;
            else %t
                obj.Sxts = Sx;
            end

            
            
            operatorMatrixIncomplete = sparse(i,j,combinedVec); %Build sparse matrix
            operatorMatrix = operatorMatrixIncomplete./dx'; %Not currently necessary, can just do /dx (scalar)

            operatorMatrix = Sx*operatorMatrix;
                
        end
        function operatorMatrix = dx2D(obj, originSubgrid)
            %Suck some things into variables for ease of use
            domain = obj.domain;  %#ok<*PROPLC,*PROP>
            BC = domain.BC;
            
            %Note: may be useful in future: for s->t matrices:
            %Ni or num T points in x direction = NxT*NyS
            
            %% TODO: handle different args for user friendliness
            
            %% Build  sparse matrix as block matrix out of 1D matrices
            operatorMatrix1D = obj.dx1D(originSubgrid);
            
            %Number of "blocks" = number of rows of in computational domain
            %(NxS(2) or NyS)
            
            NyS = domain.NxS(2);
            
            repeatedMatrix = repmat(operatorMatrix1D,1,NyS);
            repeatedMatrixCell = mat2cell(repeatedMatrix, size(operatorMatrix1D,1), repmat(size(operatorMatrix1D,2),1,NyS));
            
%             %Build argument list for blkdiag
%             diagArgs = cell(1,NyS);
%             for blockNum = 1:NyS
%                diagArgs{blockNum} = operatorMatrix1D; 
%             end
%             
            
            
            %Left and right BC already completed by 1D dx
            
%             %Take care of top and bottom boundary conditions
%             %(Note that symS and antiSymT, aka PMC on circle and triangle
%             %have no effect on top/bottom BCs)
%             
%             %Bottom:
%             if BC(2,1) == BCs.symT || BC(2,1) == BCs.antiSymS
% %                diagArgs{1} = zeros(size(operatorMatrix1D));
%                diagArgs = diagArgs(2:end);
%             end
%             %Top:
%             if BC(2,2) == BCs.symT || BC(2,2) == BCs.antiSymS
% %                 diagArgs{end} = zeros(size(operatorMatrix1D));
%                 diagArgs = diagArgs(1:end-1);
%             end
            
            %Complete "starter" matrix, and clean up
            operatorMatrix = blkdiag(repeatedMatrixCell{:});
        end
        
        function operatorMatrix = dy(obj,originSubgrid)
            %METHOD1 Returns d/dx of current operator matrix
            %   Detailed explanation goes here
            
            
            %Suck some things into variables for ease of use
            domain = obj.domain;  %#ok<*PROPLC,*PROP>
            BC = domain.BC;
            y = domain.y;
%% Todo
            
            if domain.numDimensions <=1
                error('Domain does not specify this many dimensions');
            end
            
            if nargin == 1 %originSubgrid arg is blank
                if ~(strcmp(obj.dyOriginSubgrid,'s') || strcmp(obj.dyOriginSubgrid,'t'))
                    error('dyOriginSubgrid has not yet been set to s or t.');
                end
                
                %swap the originSubgrid
                if obj.dyOriginSubgrid == 's'
                    obj.dyOriginSubgrid = 't'; 
                elseif obj.dyOriginSubgrid == 't'
                    obj.dyOriginSubgrid = 's'; 
                end
            elseif originSubgrid == 0
                if ~(strcmp(obj.dyOriginSubgrid,'s') || strcmp(obj.dyOriginSubgrid,'t'))
                    error('originSubgrid has not yet been set to s or t.');
                end
                %do not swap
            elseif strcmpi(originSubgrid, 's') || strcmpi(originSubgrid, 'e')
                obj.dyOriginSubgrid = 's'; 
                if ~isempty(obj.dyst)
                    operatorMatrix = obj.dyst;
                    return;
                end
                %otherwise, make the matrix and return it
             elseif strcmpi(originSubgrid, 't') || strcmpi(originSubgrid, 'h')
                obj.dyOriginSubgrid = 't'; 
                if ~isempty(obj.dyts)
                    operatorMatrix = obj.dyts;
                    return;
                end
            end
                
            %i is first dimension, or row # (also "destination dimension" 
            %of differential matrix)
            %j is second dimension, or column # (also "origin dimension" of
            %differential matrix)
            
            if strcmp(originSubgrid,'s')
                s = true;
            elseif strcmp(originSubgrid,'t')
                s = false;
            end
            t = ~s;
             
            
            %% The meat
            NxT = domain.NxT(1);
            NyT = domain.NxT(2);
            NxS = domain.NxS(1);
            NyS = domain.NxS(2);
            
            
            if s
                Ni = NxS*NyT; %t
                Nj = NxS*NyS; %s
            else %t
                Ni = NxS*NyS; %s
                Nj = NxS*NyT; %t
            end
           
       
            %Set up indices for sparse matrix
            %(nothing intuitive, extracted from writing matrices by hand)
            if ((BC(2,1) == BCs.antiSymS || BC(2,1) == BCs.symT) && s) ||...
                ((BC(2,1) == BCs.antiSymT || BC(2,1) == BCs.symS) && t)%PEC, Ones on main diagonal
                i11 = 1;
                j11 = 1;
                i21 = 1 + NxS;
                j21 = 1;
            elseif ((BC(2,1) == BCs.symS || BC(2,1) == BCs.antiSymT) && s) ||...
                    ((BC(2,1) == BCs.symT || BC(2,1) == BCs.antiSymS) && t)%Ones on upper diagonal
                i11 = 1;
                j11 = 1+ NxS;
                i21 = 1;
                j21 = 1;
            end
            
            if ((BC(2,2) == BCs.antiSymS || BC(2,2) == BCs.symT) && s) || ...
                    ((BC(2,2) == BCs.antiSymT || BC(2,2) == BCs.symS) && t)%PEC, 
                i1end = Ni - NxS;
                j1end = Nj;
                i2end = Ni;
                j2end = Nj;
            elseif ((BC(2,2) == BCs.symS || BC(2,2) == BCs.antiSymT) && s)||...
                    ((BC(2,2) == BCs.symT || BC(2,2) == BCs.antiSymS) && t)%PMC
                i1end = Ni;
                j1end = Nj;
                i2end = Ni;
                j2end = Nj - NxS; 
            end
            
            i1 = i11:i1end;
            j1 = j11:j1end;
            i2 = i21:i2end;
            j2 = j21:j2end;

            %Build the first diagonal for this matrix
            vec1 = ones(size(i1)); 
            %Build the second diagonal for this matrix
            vec2 = -1*ones(size(i2));  
            
%             %% Cases for left boundary (editing diagonals):
% 
%             if BC(1,1) == BCs.antiSymS || BC(1,1) == BCs.symT %Equivalent to PEC on left boundary
%                 vec1(1:NxS:end) = 0;
%                 vec2(1:NxS:end) = 0;        
%             end
%             
%             %else no change
% 
%             %% Cases for right boundary (editing diagonals):
% 
%             if BC(1,2) == BCs.antiSymS || BC(1,2) == BCs.symT %Equivalent to PEC on right boundary
%                 vec1(NxS:NxS:end) = 0;
%                 vec2(NxS:NxS:end) = 0;   
%             end
            
            %else no change
            
            %% Cases for bottom boundary (editing diagaonals):
            
            if BC(2,1) == BCs.symS 
                if s
                    vec2(1:NxS) = sqrt(2)*vec2(1:NxS); %Multiplication necessary, in case some are already 0 (but most should be -1)
                else
                    vec1(1:NxS) = sqrt(2)*vec1(1:NxS);
                end
            elseif BC(2,1) == BCs.symT
                if s
                    vec1(1:NxS) = sqrt(2)*vec1(1:NxS); %Multiplication necessary, in case some are already 0 (but most should be 1)
                else
                    vec2(1:NxS) = sqrt(2)*vec2(1:NxS);
                end
%             elseif BC(2,1) == BCs.antiSymS
%                 vec2(1:NxS) = 0;
%             elseif BC(2,1) == BCs.antiSymT
%                 vec1(1:NxS) = 0;
            end
            
            %% Cases for top boundary (editing diagonals):
            
            if BC(2,2) == BCs.symS 
                if s
                    vec1(end-(NxS-1):end) = sqrt(2)*vec1(end-(NxS-1):end); 
                else
                    vec2(end-(NxS-1):end) = sqrt(2)*vec2(end-(NxS-1):end); 
                end
                %Multiplication necessary, in case some are already 0 (but most should be 1)
            elseif BC(2,2) == BCs.symT
                if s
                    vec2(end-(NxS-1):end) = sqrt(2)*vec2(end-(NxS-1):end); 
                else
                    vec1(end-(NxS-1):end) = sqrt(2)*vec1(end-(NxS-1):end); 
                end
                %Multiplication necessary, in case some are already 0 (but most should be -1)
%             elseif BC(2,2) == BCs.antiSymS
%                 vec1(end-(NxS-1):end) = 0;
%             elseif BC(2,2) == BCs.antiSymT
%                 vec2(end-(NxS-1):end) = 0;
            end
            
            %% Combine left and right diagonals into one matrix

            i = [i1,i2];
            j = [j1,j2];
            combinedVec = [vec1,vec2]'; %Create one long vector
            
            %% But first! Modify matrix based on local dx

             %Index of first s point in x vector
            if BC(2,1) == BCs.antiSymT || BC(2,1) == BCs.symS %Start on s point in x direction
                s1 = 1; %Not used, we use -op'
                t1 = 2;
            else %Start on t point in x direction
                s1 = 2; %Not used, we use -op'
                t1 = 1;
            end

            %Build dx vector
            %This is a little difficult to read, but was necessary in
            %order to use no for-loops

            %In 1D, Nj was the number of origin subgrid points in the x
            %direction, and Ni was the number of target subgrid points in
            %the x direction
            %In 2D, Nj -> NOrigRows should now be the number of origin subgrid ROWS, and
            %Ni -> NTargRows is now the number of target subgrid ROWS
            if s
                NOrigRows = NyS;
                NOrigColumns = NxS;
                %NTargRows = NyT;
            else
                NOrigRows = NyT;
                NOrigColumns = NxS; %For our purposes; since we're doing y here
                %NTargRows = NyS;
            end
            
            if s
                y1 = y(s1:2:(NOrigRows-1)*2); 
                y2 = y(s1+2:2:(NOrigRows)*2); %Note the end index here may overshoot by 1, but matlab won't go the extra half cell because of the :2:


                if BC(2,1) == BCs.antiSymT || BC(2,1) == BCs.symS %End on an s point
                    dy1 = [];
                else %End on a t point
                    dy1 = (y(2) - y(1))*2;
                end


                if BC(2,2) == BCs.antiSymT || BC(2,2) == BCs.symS %End on an s point
                    dyend = [];
                else %End on a t point
                    dyend = (y(end) - y(end-1))*2; 
                end

            else
                y1 = y(t1:2:(NOrigRows-1)*2); 
                y2 = y(t1+2:2:(NOrigRows)*2); %Note the end index here may overshoot by 1, but matlab won't go the extra half cell because of the :2:


                if BC(2,1) == BCs.antiSymT || BC(2,1) == BCs.symS %End on an s point
                    dy1 = (y(2) - y(1))*2;
                else %End on a t point
                    dy1 = [];
                end

                if BC(2,2) == BCs.antiSymT || BC(2,2) == BCs.symS %End on an s point
                    dyend = (y(end) - y(end-1))*2; 
                else %End on an t point
                    dyend = [];
                end
            end


            dy = y2-y1; %contains the dx from surrounding originSubgrid points, around target points
            dy = [dy1 dy dyend];
            
            %1/dy2D, then rep and multiply
            dy = 1./dy;
            dy2D = repmat(dy,NOrigColumns,1) ;
            dy2D = dy2D(:);
            
            
            power = 3; %Later, make this a param
            Sy = ones (size(domain.y));
                
            if domain.PML(2,1) >0 %Bottom boundary
                imComponent = (domain.PML(2,1)*domain.dx(2)):-domain.dx(2)/2:domain.dx(2)/2;
                A = domain.PMLMag /((domain.PML(2,1)*domain.dx(2))^power);
                imComponent = 1i*(A*imComponent.^power);
                %imComponent = -1* imComponent;
                Sy(1:domain.PML(2,1)*2) = Sy(1:domain.PML(2,1)*2) + imComponent;
            end
            
            if domain.PML(2,2) >0
                imComponent = domain.dx(2)/2:domain.dx(2)/2:domain.PML(2,2)*domain.dx(1);
                A = domain.PMLMag /((domain.PML(2,2)*domain.dx(2))^power);
                imComponent = 1i*(A*imComponent.^power);
                %imComponent = -1* imComponent;
                Sy(end-domain.PML(2,2)*2+1:end) = Sy(end-domain.PML(2,2)*2+1:end) + imComponent;
            end
            
            if s
               Sy = domain.GetSubvector('t',2,Sy); 
            else
                Sy = domain.GetSubvector('s',2,Sy);
            end
            
            Sy = repmat(Sy,NOrigColumns,1) ;
            if s 
                obj.Syst = Sy;
            else
                obj.Syts = Sy;
            end
            Sy = Sy(:);
            Sy = diag(sparse(1./Sy));
            
            
            
            %% Finally, build sparse matrix
            if length(i) ~= length(j) || length(i) ~= length(combinedVec)
               disp('uh oh'); 
            end
            operatorMatrixIncomplete = sparse(i,j,combinedVec); %Build sparse matrix
            
%             tic
%             operatorMatrix = dy2D.*operatorMatrixIncomplete;
%             toc
%             tic
            operatorMatrix = bsxfun(@times, dy2D, operatorMatrixIncomplete); 
%             toc
            %Try inverting dy2D first, then multiply
            %^This final step is what takes most of the time... can we
            %speed it up?
            
            operatorMatrix = Sy*operatorMatrix;
            
            
         end
        function operatorMatrix = dz(obj,originSubgrid)
            %METHOD1 Returns d/dx of current operator matrix
            %   Detailed explanation goes here
            
            
            %Suck some things into variables for ease of use
            domain = obj.domain;  %#ok<*PROPLC,*PROP>
            BC = domain.BC;
            operatorMatrix = -1;
            if domain.numDimensions <=2
                error('Domain does not specify this many dimensions');
            end
        end
    end
end

