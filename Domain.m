classdef Domain
    %DOMAIN Collection of parameters describing computational domain
    %Padding of endpoints automatically to ensure correct boundary conditions
    
    properties
        dx 
        NxS 
        NxT 
        x %2x grid, with 0 always starting on original s point, and padded around this
        y %Cannot contain these in one matrix, as they may have different dim
        z %Could use cell, but clunky syntax
        BC 
        PML
        PMLMag
        numDimensions
        eps %defined on the 2x grid, sparse diagonal matrices
        mu  %defined on the 2x grid, sparse diagonal matrices
    end
    
    methods
        function obj = Domain(dx,Nx,BC,PML, PMLMag, numDimensions)
            %DOMAIN Construct an instance of this class
            
            obj.dx = dx;    %Vector: [dx;dy;dz]
            obj.NxS = round(Nx);   %Vector: [NxS; NyS; NzS] %Depending on BC, these may actually be up to Nx+1
            obj.NxT = round(Nx);   %Vector: [NxT; NyT; NzT]
            obj.BC = BC;    %Matrix: %x-, x+; y-, y+; z-, z+ in that order. 
                            %E.g.: BC = [BCs.antiSymS, BCs.antiSymT; BCs.symS,BCs.symT; etc.,etc.]
            obj.PML = round(PML*2)/2;  %Matrix containing the number of desired PML (full) cells on each boundary
                            %This is overlaid onto the existing domain: [-x,x;-y,y;-z,z]
                            %e.g. [20,20;0,0] for 2D, with no PML on top and bottom
                            %Can handle half cells; use multiples of 0.5
            obj.PMLMag = PMLMag; %A scalar designating the max value of the imaginary component 
                                 %at the edge of all boundaries
            
            obj.numDimensions = numDimensions;
            
            if (obj.numDimensions >=1)
                
                x1 = 0;                         %Standard start point on 2x grid (at first s point)

                xN = dx(1)/2*(Nx(1)*2-1);       %Standard end point on a 2x grid (at final t half point)

                if BC(1,1) == BCs.antiSymS      %Need to add a T before the first S point, so the new first S = 0
                    x1 = -dx(1)/2;
                    obj.NxT(1) = obj.NxT(1) +1;
                    
                %elseif BC(1,1) == BCs.symS %No change, we get sqrt(2) condition, but no padding needed on left
                %elseif BC(1,1) == BCs.antiSymT %There is a T point on the left now, but it is 0, so we don't add anything
                
                elseif BC(1,1) == BCs.symT      %Add a T point on the left
                    x1 = -dx(1)/2;
                    obj.NxT(1) = obj.NxT(1) +1;
                end

                if BC(1,2) == BCs.antiSymS      %S point on end, but == 0, so don't add anything
                elseif BC(1,2) == BCs.symS      %Add S point on end
                    xN = xN + dx(1)/2;
                    obj.NxS(1) = obj.NxS(1) +1;
                elseif BC(1,2) == BCs.antiSymT  %T point on end, but == 0, so pad with extra S point
                    xN = xN + dx(1)/2;
                    obj.NxS(1) = obj.NxS(1) +1;
                    
                %elseif BC(1,2) == BCs.symT %T point on end alread, do nothing
                
                end

                obj.x = x1:dx(1)/2:xN;          %Creates x as list of positions, in 2x resolution, 
                                                %up to Nx*2+1 points in total
                if mod(PML(1,1),0.5)~=0 || mod(PML(1,2),0.5)~=0
                    error('Invalid PML, must be positive multiple of 0.5');
                end
                if length(obj.x) < obj.PML(1,1)*2 + obj.PML(1,2)*2
                    error('Desired PMLs would take up entire domain (x)!')
                    
                end
                
          %This may be a mathematically incorrect way of doing it...
%                 %% Add PMLs in x vector via complex coordinate stretching
%                 if obj.PML(1,1) > 0 %Left boundary 
%                     % TODO: Should increase as s(x)/omega so all frequencies will
%                     % decay at same rate
%                     % s(x) should be a cubic function for starters
%                     % Increase up to max value: PMLMag
%                     
%                     imComponent = (obj.PML(1,1)*dx(1)):-dx(1)/2:dx(1)/2;    % E.g.: 20:0.5:0.5
%                     A = obj.PMLMag /((obj.PML(1,1)*dx(1))^3);               % Calculate scaling factor to end at PMLMag 
%                     imComponent = 1i*(A*imComponent.^3);                    % Generate cubically-increasing imaginary component
%                     imComponent = -1* imComponent;
%                     obj.x(1:obj.PML(1,1)*2) = obj.x(1:obj.PML(1,1)*2) + imComponent; 
%                                                                             % ^Overlay generated complex coordinates onto
%                                                                             % left side of existing x vector
%                 end
%                 
%                 if obj.PML(1,2) > 0 %Right boundary
%                     imComponent = dx(1)/2:dx(1)/2:obj.PML(1,2)*dx(1);
%                     A = obj.PMLMag /((obj.PML(1,2)*dx(1))^3);
%                     imComponent = 1i*(A*imComponent.^3);
%                     imComponent = -1* imComponent;
%                     obj.x(end-obj.PML(1,2)*2+1:end) = obj.x(end-obj.PML(1,2)*2+1:end) + imComponent;
%                 end
                
            end
            if (obj.numDimensions >=2)
                y1 = 0;
                yN = dx(2)/2*(Nx(2)*2-1);
                
                if BC(2,1) == BCs.symT || BC(2,1) == BCs.antiSymS  %Need to add(pad) a t point to bottom
                    y1 = -dx(2)/2;
                    obj.NxT(2) = obj.NxT(2) +1;
                end
                if BC(2,2) == BCs.symS || BC(2,2) == BCs.antiSymT  %Need to add(pad) an s point to top
                    yN = yN + dx(2)/2;
                    obj.NxS(2) = obj.NxS(2) +1;
                end
                obj.y = y1:dx(2)/2:yN;
                
            %This may be a mathematically incorrect way of doing it...
%                 %% Add PMLs in y vector via complex coordinate stretching
%                 if mod(PML(2,1),0.5)~=0 || mod(PML(2,2),0.5)~=0
%                     error('Invalid PML, must be positive multiple of 0.5');
%                 end
%                 if length(obj.y) < obj.PML(2,1)*2 + obj.PML(2,2)*2
%                     error('Desired PMLs would take up entire domain (y)!')
%                     
%                 end
%                 
%                 if obj.PML(2,1) >0
%                     imComponent = (obj.PML(2,1)*dx(2)):-dx(2)/2:dx(2)/2;
%                     A = obj.PMLMag /((obj.PML(2,1)*dx(2))^3);
%                     imComponent = 1i*(A*imComponent.^3);
%                     imComponent = -1* imComponent;
%                     obj.y(1:obj.PML(2,1)*2) = obj.y(1:obj.PML(2,1)*2) + imComponent;
%                 end
%                 if obj.PML(2,2) >0
%                     imComponent = dx(2)/2:dx(2)/2:obj.PML(2,2)*dx(1);
%                     A = obj.PMLMag /((obj.PML(2,2)*dx(2))^3);
%                     imComponent = 1i*(A*imComponent.^3);
%                     imComponent = -1* imComponent;
%                     obj.y(end-obj.PML(2,2)*2+1:end) = obj.y(end-obj.PML(2,2)*2+1:end) + imComponent;
%                 end
                
            end
            if (obj.numDimensions >=3)
                error('3d not implemented')
%                 z1 = 0;
%                 zN = dx(3)/2*(Nx(3)*2-1);
            end
            obj = SetConstantPropertyGrid(obj, 'eps', 1);
            obj = SetConstantPropertyGrid(obj, 'mu', 1);
        end
        
        function grid = GetXGrid(obj, dimension)
            if obj.numDimensions == 1
                if dimension == 1
                    grid = obj.x;
                else
                    error('Dimension not supported');
                end
            elseif obj.numDimensions == 2
                if dimension == 1
                    [grid,~] = meshgrid(obj.x, obj.y);
                elseif dimension == 2
                    [~,grid] = meshgrid(obj.x, obj.y);
                else
                    error('Dimension not supported');
                end
            else
                error('Dimension not supported');
            end
        end
        
        function obj = SetConstantPropertyGrid(obj, property, value)
            %Accepts property ('eps' or 'mu')
            %Creates a 2x grid 
            %Sets each grid point to 'value'
            %Saves grid to obj
            %This should output an x vector (in 1D or 2D)
            %Input corresponds to the x y coordinate vector [x;y]
            
            propMatrix = value*ones(size(obj.GetXGrid(obj.numDimensions)));
%             propMatrix = propMatrix(:);
%             propMatrix = diag(sparse(propMatrix));

            if strcmp(property, 'eps')
               obj.eps = propMatrix; 
            elseif strcmp(property, 'mu')
               obj.mu = propMatrix;
            else
                error('Properties supported: mu and eps');
            end
        end
        
        function propertySubgrid = GetPropertySubgrid(obj, property, diagonalize)
            %Accepts property ('eps' or 'mu')
            %Accepts subgrid ('x', 'y', or 'z')
            %Returns grid of, e.g. eps_x, or mu_z
%             if all(obj.BC == BCs.antiSymS)
                %epszz: 1:2:end
                %
            %if all are antisymmetric t:
                %epszz: 1:2:end  
                
            Nx2 = obj.NxS(1) + obj.NxT(1);
            Ny2 = obj.NxS(2) + obj.NxT(2);
            %Ez mode
            if strcmp(property,'eps_zz')
                propertySubgrid = obj.eps(1:2:end,1:2:end);
            elseif strcmp(property,'mu_xx')
                propertySubgrid = obj.mu(1:2:Nx2,2:2:Ny2);
            elseif strcmp(property,'mu_yy')
                propertySubgrid = obj.mu(2:2:Nx2,1:2:Ny2);
                
            %Hz mode
            elseif strcmp(property,'eps_xx')     
                propertySubgrid = obj.eps(2:2:Nx2,1:2:Ny2);
            elseif strcmp(property,'eps_yy') 
                propertySubgrid = obj.eps(1:2:Nx2,2:2:Ny2);
            elseif strcmp(property,'mu_zz') 
                propertySubgrid = obj.mu(2:2:Nx2,2:2:Ny2);
            end
            
            if diagonalize
                propertySubgrid = propertySubgrid';
                propertySubgrid = propertySubgrid(:);
                propertySubgrid = diag(sparse(propertySubgrid));
            end
            %Subgrid is not necessarily x or y or s or t... it is xx, yy,
            %zz, and this is different based on whether mu or eps is
            %requested - see rumpf implementation/ slides for a reminder
        end
        function PlotGrid(obj,property,plotGrid)
            figure;
            Nx2 = obj.NxS(1) + obj.NxT(1);
            Ny2 = obj.NxS(2) + obj.NxT(2);
            if strcmp(property,'eps')
                imagesc(obj.eps);
                axis xy;
                grid off;
                axis equal;
                axis image;
                title('eps (2x resolution)');
                colorbar; colormap(gray);
                if plotGrid
                     xticks(1:Nx2);
                     yticks(1:Ny2);

                     hold on
                    g_y=0.5:2:Ny2; % user defined grid Y [start:spaces:end]
                    g_x=0.5:2:Nx2; % user defined grid X [start:spaces:end]
                    for i=1:length(g_x)
                       plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'k:', 'Color', 'red') %y grid lines
                       hold on    
                    end
                    for i=1:length(g_y)
                       plot([g_x(1) g_x(end)],[g_y(i) g_y(i)],'k:','Color', 'red') %x grid lines
                       hold on    
                    end 
                end
            elseif strcmp(property,'mu')
                imagesc(obj.mu);
                axis xy;
                grid off;
                axis equal;
                axis image;
                title('mu (2x resolution)');
                colorbar; colormap(gray);
                if plotGrid
                     xticks(1:Nx2);
                     yticks(1:Ny2);

                     hold on
                    g_y=0.5:2:Ny2; % user defined grid Y [start:spaces:end]
                    g_x=0.5:2:Nx2; % user defined grid X [start:spaces:end]
                    for i=1:length(g_x)
                       plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'k:', 'Color', 'red') %y grid lines
                       hold on    
                    end
                    for i=1:length(g_y)
                       plot([g_x(1) g_x(end)],[g_y(i) g_y(i)],'k:','Color', 'red') %x grid lines
                       hold on    
                    end 
                end
            end
        end
        function obj = AddRectangle(obj, property, value, boundaries)
            %Adds to 2x grid!
            %Accepts left, right, top, and bottom boundaries ([l,r;t,b])
            %in terms of x and y coordinate vectors
            %Overwrites the positions within that rectangle (inclusive) to a
            %the value specified
            l = boundaries(1,1);
            r = boundaries(1,2);
            t = boundaries(2,1);
            b = boundaries(2,2);
            if t>b || l>r
                error('Bad boundaries')
            end
            
            Nx2 = obj.NxS(1) + obj.NxT(1);
            Ny2 = obj.NxS(2) + obj.NxT(2);
            
            if strcmp(property,'eps')
                obj.eps(t:b,l:r) = value;
            elseif strcmp(property,'mu')
                obj.mu(t:b,l:r) = value;
            
            %Ez mode
            elseif strcmp(property,'eps_zz')
                %obj.eps(1:2:Nx2,1:2:Ny2) = ;
            elseif strcmp(property,'mu_xx')
                %propertySubgrid = obj.mu(1:2:Nx2,2:2:Ny2);
            elseif strcmp(property,'mu_yy')
                %propertySubgrid = obj.mu(2:2:Nx2,1:2:Ny2);
                
            %Hz mode
            elseif strcmp(property,'eps_xx')     
                %propertySubgrid = obj.eps(2:2:Nx2,1:2:Ny2);
            elseif strcmp(property,'eps_yy') 
                %propertySubgrid = obj.eps(1:2:Nx2,2:2:Ny2);
            elseif strcmp(property,'mu_zz') 
                %propertySubgrid = obj.mu(2:2:Nx2,2:2:Ny2);
            end
           
            
        end
        
        
        function obj = AddRectangleInterpolatedEdges(obj,property,value,boundaries)
            %Adds to 2x grid!
            %Accepts left, right, top, and bottom boundaries ([l,r;t,b])
            %in terms of x and y coordinate vectors
            %Overwrites the positions within that rectangle (inclusive) to a
            %the value specified
            l = boundaries(1,1);
            r = boundaries(1,2);
            t = boundaries(2,1);
            b = boundaries(2,2);
            if t>b || l>r
                error('Bad boundaries')
            end
            
            fd = @(x) 1./(1 + exp(x));
            ws = obj.dx(1)/4;  % Smoothing distance
%             x1 = position(1);
%             y1 = position(2);
%             R = radius;
            X = obj.GetXGrid(1);
            Y = obj.GetXGrid(2);
            
            %obj.eps = obj.eps + (value-obj.eps)*(fd((Y-(b))/ws)-fd((Y-(t))/ws))*(fd((X-(l))/ws)-fd((X-(r))/ws));  % Add bus wg
            
            Nx2 = obj.NxS(1) + obj.NxT(1);
            Ny2 = obj.NxS(2) + obj.NxT(2);
            
            if strcmp(property,'eps')
                %obj.eps(t:b,l:r) = value;
                obj.eps = obj.eps + (value-obj.eps).*(fd((Y-b)/ws)-fd((Y-t)/ws)).*(fd((X-r)/ws)-fd((X-l)/ws));
            elseif strcmp(property,'mu')
                %obj.mu(t:b,l:r) = value;
                obj.mu = obj.mu + (value-obj.mu).*(fd((Y-b)/ws)-fd((Y-t)/ws)).*(fd((X-r)/ws)-fd((X-l)/ws));
            
            %Ez mode
            elseif strcmp(property,'eps_zz')
                %obj.eps(1:2:Nx2,1:2:Ny2) = ;
            elseif strcmp(property,'mu_xx')
                %propertySubgrid = obj.mu(1:2:Nx2,2:2:Ny2);
            elseif strcmp(property,'mu_yy')
                %propertySubgrid = obj.mu(2:2:Nx2,1:2:Ny2);
                
            %Hz mode
            elseif strcmp(property,'eps_xx')     
                %propertySubgrid = obj.eps(2:2:Nx2,1:2:Ny2);
            elseif strcmp(property,'eps_yy') 
                %propertySubgrid = obj.eps(1:2:Nx2,2:2:Ny2);
            elseif strcmp(property,'mu_zz') 
                %propertySubgrid = obj.mu(2:2:Nx2,2:2:Ny2);
            end
            
            
            
        end

        
        function obj = AddDiskInterpolatedEdges(obj, property, value, position, radius)
            fd = @(x) 1./(1 + exp(x));
            ws = obj.dx(1)/4;  % Smoothing distance
            x1 = position(1);
            y1 = position(2);
            R = radius;
            X = obj.GetXGrid(1);
            Y = obj.GetXGrid(2);
            if strcmp(property,'eps')
                obj.eps = obj.eps +(value-obj.eps).*fd((sqrt((X-x1).^2 + (Y-y1).^2) - R)/ws); % Draw disk with outer radius - ring 1;
            elseif strcmp(property,'mu')
                obj.mu = obj.mu +(value-obj.mu).*fd((sqrt((X-x1).^2 + (Y-y1).^2) - R)/ws); % Draw disk with outer radius - ring 1;
            end
        end
       
        
        function obj = AddRing(obj, property, value, position, radius,width)
            obj = obj.AddDiskInterpolatedEdges(property,value,position,radius+width/2);
            obj = obj.AddDiskInterpolatedEdges(property,1,position,radius-width/2);
        end
        
        %Generic: not for x/y in particular (pass in the desired full
        %vector), note fullVector must be the corresponding vector!
               
        %Supports s, t, u subgrids
        function subvec = GetSubvector(obj, originSubgrid, dimension, fullVector)
            [s,t,u]= obj.ParseSubgridString(originSubgrid);
            
            if dimension == 1
                if obj.numDimensions < 1
                    error('How did you manage to make a 0D Domain?')
                end
                
                if s || u
                    if obj.BC(1,1) == BCs.antiSymS || obj.BC(1,1) == BCs.symT     
                        i1 = 2;
                    else
                        i1 = 1;
                    end

                    subvec = fullVector(i1:2:obj.NxS(1)*2);

                elseif t
                    if obj.BC(1,1) == BCs.antiSymS || obj.BC(1,1) == BCs.symT     
                        i1 = 1;
                    else
                        i1 = 2;
                    end

                    subvec = fullVector(i1:2:obj.NxT(1)*2);
                end
            elseif dimension == 2
                if obj.numDimensions < 2
                    error('Entered value for dimension not supported');
                end
           
                if s
                    if obj.BC(2,1) == BCs.symT || obj.BC(2,1) == BCs.antiSymS  %Need to add(pad) a t point to bottom
                        i1 = 2;
                    else
                        i1 = 1;
                    end

                    subvec = fullVector(i1:2:obj.NxS(2)*2);
                        
                elseif t || u
                    if obj.BC(2,1) == BCs.symT || obj.BC(2,1) == BCs.antiSymS  %Need to add(pad) a t point to bottom
                        i1 = 1;
                    else
                        i1 = 2;
                    end

                    subvec = fullVector(i1:2:obj.NxT(2)*2);

                end
            end
        end
            
        function x_subvec = GetXSubvector(obj, originSubgrid, dimension)
            if dimension == 1
                x_subvec = GetSubvector(obj, originSubgrid, dimension, obj.x);
            elseif dimension == 2
                x_subvec = GetSubvector(obj, originSubgrid, dimension, obj.y);
            end
        end
        
        %?
        function subgrid = GetSubgrid(obj, originSubgrid, dimension, fullGrid)
            %This may be trickier than it first seems!
            if obj.numDimensions == 1
                if dimension == 1
                    subgrid = GetSubvector(obj, originSubgrid, dimension, obj.x);
                else
                    error ('dimension not supported');
                end
            elseif obj.numDimensions == 2
                if dimension == 1
                    
                elseif dimension == 2
                    
                else
                   error('dimension not supported') 
                end
                
            else
                error('dimension not supported')
            end
        end
        
        %Supports s,t,u
        function x_subgrid = GetXSubgrid(obj,originSubgrid, dimension)
           x_subgrid = GetXSubgridUnwrapped(obj, originSubgrid, dimension);
           
           if obj.numDimensions == 1
               return;
           end
           %Now reshape it!
           [s,t,u] = obj.ParseSubgridString(originSubgrid);
           if s
               x_subgrid = reshape(x_subgrid, obj.NxS(1), obj.NxS(2));
           elseif t
               x_subgrid = reshape(x_subgrid, obj.NxT(1), obj.NxT(2));
           elseif u
               x_subgrid = reshape(x_subgrid, obj.NxS(1), obj.NxT(2));
           end
        end
        
        function xIndex = GetXIndex(obj, xPosition, dimension)
            if dimension ==1
                [ ~, xIndex ] = min( abs( obj.x-xPosition ) );
            elseif dimension == 2
                [ ~, xIndex ] = min( abs( obj.y-xPosition ) );
            else
                error('Dimension not supported')
            end
        end
        
        
        function [s,t,u] = ParseSubgridString(obj, subgrid)
            
            s = false;
            t = false;
            u = false;
           
            if strcmpi(subgrid, 's') || strcmpi(subgrid, 'e')
                s = true; 
            elseif strcmpi(subgrid, 't') || strcmpi(subgrid, 'h')
                t = true;
            elseif strcmpi(subgrid, 'u')
                u = true;
            else 
               error('Entered value (%s) for originSubgrid argument is not supported', subgrid); 
            end
        end
        
        
        
        
        %Note a Mx gaussian source will be defined on the u points, and an
        %My gaussian source will be defined on the t points
        %This will be the same as:
        %Mx: GetXSubgridUnwrapped(obj, 'u', 1) 
        %My: GetXSubgridUnwrapped(obj, 't', 2)
        
        %This now will support s, t, and u
        function x_subgrid = GetXSubgridUnwrapped(obj, originSubgrid, dimension)
            
            x_subvec = obj.GetXSubvector(originSubgrid, dimension);
            
            [s,t,u] = obj.ParseSubgridString(originSubgrid);
            
            if dimension == 1
                if obj.numDimensions == 1
                    x_subgrid = x_subvec;
                elseif obj.numDimensions == 2
                    %We need to replicate this x vector by the number of rows.
                    if s || t
                        numRows = obj.NxS(2); %Number of s/t rows!
                    elseif u
                        numRows = obj.NxT(2); %Number of u rows!
                    end
                        
                    x_subgrid = repmat(x_subvec, 1, numRows);
                else
                    error('Num dimensions not supported');
                end
            elseif dimension == 2
                if obj.numDimensions ~= 2
                    error('Num dimensions not supported');
                end
                %Need to provide a yvec, but with y(1) repeated once
                %for each column
                if s || u
                    numColumns = obj.NxS(1);
                elseif t
                    numColumns = obj.NxT(1); 
                end
                x_subgrid = repelem(x_subvec, numColumns);
            end
        end
        
        %Supports s,t,u
        function indexVec = GetWrappedIndex(obj, x, originSubgrid)
            if obj.numDimensions == 1
                if ~all(size(x) == [1 1])
                    error('1D domain requires 1D x position input for delta vector');
                end
                
                xsubvec = obj.GetXSubvector(originSubgrid, 1);
                %Find the nearest value in x vector to inputted x
                [ ~, indexVec(1) ] = min( abs( xsubvec-x ) );
                %Above, ~ is the dx between x and the closest val in obj.x
                
            elseif obj.numDimensions == 2
                if ~all(size(x) == [2 1])
                    error('2D domain requires 2D x position input for delta vector (did you forget the semicolon?)');
                end
                
                %Find the nearest value in x vector to inputted x
                xsubvec = obj.GetXSubvector(originSubgrid,1);
                ysubvec = obj.GetXSubvector(originSubgrid,2);
                
                [ ~, indexVec(1) ] = min( abs( xsubvec-x(1) ) );
                %Above, ~ is the dx between x and the closest val in obj.x
                
                %Find the nearest value in y vector to inputted y (x(2))
                [ ~, indexVec(2) ] = min( abs( ysubvec-x(2) ) );
                
            else
                error('Num dimensions not supported');
            end
        end
        
        %Only works for s currently
        function totalIndex = GetUnwrappedIndex(obj, x, originSubgrid)
            indexVec = GetWrappedIndex(obj, x, originSubgrid);
            
            if obj.numDimensions == 1
                totalIndex = indexVec;
                
            elseif obj.numDimensions == 2
                rowLength = obj.NxS(1);
                totalIndex = indexVec(1) + rowLength*(indexVec(2)-1);
            else
                error('Num dimensions not supported');
            end
        end
        
        %TODO: make it work for both
        %Currently only functioning for 's' subgrid!
        function delta = Delta(obj, x,originSubgrid)
            %This should output an x vector (in 1D or 2D)
            %Input corresponds to the x y coordinate vector [x;y]
            
            if obj.numDimensions == 1
                if ~all(size(x) == [1 1])
                    error('1D domain requires 1D x position input for delta vector');
                end

                xIndex = obj.GetUnwrappedIndex(x,originSubgrid);
                    
                delta = zeros(size(obj.x));
                delta(xIndex) = 1;
                
            elseif obj.numDimensions == 2
                if ~all(size(x) == [2 1])
                    error('2D domain requires 2D x position input for delta vector (did you forget the semicolon?)');
                end
                
                totIndex = obj.GetUnwrappedIndex(x,originSubgrid);
                
                delta = zeros(size(GetXSubgridUnwrapped(obj,originSubgrid,1)));
                delta(totIndex) = 1;
                
            else
                error('Num dimensions not supported');
            end
            delta = delta(:);
        end
        
        %Currently only functioning for 's' subgrid
        %x: coordinates of center of Gaussian
        %dimension: direction of Gaussian (Note this is defined only at one
        %"slice" in either the x or y direction)
        %originSubgrid: the grid to define the gaussian on
        %a,c: shape parameters
        
        function gaussian = Gaussian(obj, x, originSubgrid, dimension, a, c)
            %Defined as a*exp(-(x-x_param)^2/(2c^2))
            if obj.numDimensions == 1
                if ~all(size(x) == [1 1])
                    error('1D domain requires 1D x position input for delta vector (did you forget to clear?)');
                end
                if dimension > 1 || dimension < 1
                    error('Dimension not supported')
                end
                
                %xIndex = obj.GetUnwrappedIndex(x);
                    
%               gaussian = zeros(size(obj.x));
                xSubgrid = obj.GetXSubvector(originSubgrid,1);
                gaussian = a*exp(-((xSubgrid-x).^2)/(2*c^2));
                
            elseif obj.numDimensions == 2
                if ~all(size(x) == [2 1])
                    error('2D domain requires 2D x position input for delta vector (did you forget the semicolon?)');
                end
                
                if dimension == 1
                    
                    gaussianVec = obj.GetXSubvector(originSubgrid, 1);
                    gaussianVec = a*exp(-((gaussianVec - x(1)).^2)/(2*c^2));
                    
                    gaussian = zeros(size(GetXSubgrid(obj,originSubgrid,1)));
                    indexVec = GetWrappedIndex(obj, x, originSubgrid);
                    
                    gaussian(:,indexVec(2)) = gaussianVec;
                    
                elseif dimension == 2 
                    gaussianVec = obj.GetXSubvector(originSubgrid, 2);
                    gaussianVec = a*exp(-((gaussianVec - x(2)).^2)/(2*c^2));
                    
                    gaussian = zeros(size(GetXSubgrid(obj,originSubgrid,2)));
                    indexVec = GetWrappedIndex(obj, x, originSubgrid);
                    
                    gaussian(indexVec(1),:) = gaussianVec;
                else
                    error('Dimension not supported')
                end
                
                
                
                
            else
                error('Num dimensions not supported');
            end
            %Unwrap
            gaussian = gaussian(:);
        end
        
        function midpoint = GetMidpoint(obj)
           
            midpoint(1) = obj.x(round(length(obj.x)/2));
            if obj.numDimensions > 1
                syvec = obj.GetXSubvector('s', 2);
                midpoint(2) = syvec(round(length(syvec)/2));; 
            end
        end
        
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

