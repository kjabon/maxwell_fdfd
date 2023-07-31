function scalarModes2D(mu, eps, dx, NxIn, boundaryConditions, pmls, pmlMag)
    %demo3
    %2D modes of cavity, with zero boundaries all around
    clc;

    %mu and eps need to be sparse diagonal matrices
    
%     mu = 1;
%     eps = 1; %Later these need to be 1d arrays, i.e. an expanded ND matrix
%     dx = [1, 1];
%     NxIn = [300, 300];
%     boundaryConditions = [BCs.antiSymS, BCs.antiSymS; BCs.antiSymS, BCs.antiSymS];
%     pmls = [0,0;0,0];
%     pmlMag = 0; 
    numDimensions = 2;

    %Create domain and diffOperator instances
    myDomain = Domain(dx,NxIn,boundaryConditions,pmls,pmlMag,numDimensions);
    D = DifferentialOperator(myDomain);
    dxSquared = D.dxts*D.dxst;
    dySquared = D.dyts*D.dyst;
    delSquared = dxSquared+dySquared;
    A = -1*(1/mu)*(1/eps)*delSquared;

    [V,eigenvalues] = eigs(A,5,'smallestabs');
    % plot(V(:,1));

    Nx = myDomain.NxS(1);
    Ny = myDomain.NxS(2);
    V1 = reshape(V(:,1),Nx, Nx);
    %V1 = abs(V1);
     V2 = reshape(V(:,2),Nx(1), Ny);
     V3 = reshape(V(:,3),Nx(1), Ny);
     V4 = reshape(V(:,4),Nx(1), Ny);
     V5 = reshape(V(:,5),Nx(1), Ny);
     %V6 = reshape(V(:,6),Nx(1), Ny);
    imagesc(V1)
    pbaspect([1 1 1])

    figure
    %V2 = abs(V2);
    imagesc(V2)
    pbaspect([1 1 1])

    figure
    %V3 = abs(V3);
    imagesc(V3)
    pbaspect([1 1 1])

    figure
    %V4 = abs(V4);
    imagesc(V4)
    pbaspect([1 1 1])

    figure
    % V5 = abs(V5);
    imagesc(V5)
    pbaspect([1 1 1])

    % figure
    % V6 = abs(V6);
    % imagesc(V6)
    % pbaspect([1 1 1])

    % figure;
    % Vtot = V1+V2+V3+V4+V5;%+V6;
    % imagesc(Vtot);
    % pbaspect([1 1 1]);

end