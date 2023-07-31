%Test differential operator

BCS = [BCs.antiSymS BCs.antiSymT BCs.symS BCs.symT];
numDimensions = 2;
sizes = 4;
visualization = true;
if visualization
    close all;
    colormap(jet(128));
    
end
for size1 = sizes
        disp(sprintf('Grid: %d x %d', size1,size1)); %#ok<DSPS>
        for BC1 = BCS
            
            for BC2 = BCS
                if numDimensions == 2
                    for BC3 = BCS
                        for BC4 = BCS
                            disp([BC1 BC2;BC3 BC4]);
                            myDomain = Domain([1;1],[size1;size1],[BC1, BC2; BC3, BC4],[2,2;2,2],2,2)
                            disp('dx:')
                            disp(myDomain.dx)
                            disp('NxS:')
                            disp(myDomain.NxS)
                            disp('NxT:')
                            disp(myDomain.NxT)
                            disp('Boundary conditions:')
                            disp(myDomain.BC)
                            disp('Num PML cells:')
                            disp(myDomain.PML)
                            
                            disp('x vector:')
                            disp(myDomain.x)
                            disp('y vector:')
                            disp(myDomain.y)
                            
                            
                            D = DifferentialOperator(myDomain)
                            if BC1 == BCs.antiSymS && BC2 == BCs.antiSymT && BC3 == BCs.symS && BC4 == BCs.symT
                                disp('break')
                            end
                            if visualization
                                if size1 <= 10
                                    disp('d/dx matrix from s to t');
                                    spy(D.dxst)
                                    set(gcf,'color','w');
                                    disp('d/dx matrix from t to s');
                                    spy(D.dxts)
                                    set(gcf,'color','w');
                                    disp('d/dy matrix from s to t');
                                    spy(D.dyst)
                                    set(gcf,'color','w');
                                    disp('d/dy matrix from t to s');
                                    spy(D.dyts)
                                    set(gcf,'color','w');
%                                     imagesc(D.dxst);
%                                     pause(0.1);
%                                     imagesc(D.dyst);
%                                     pause(0.1);
                                elseif size1 <= 16
                                    imagesc(D.dyst);
                                    pause(0.05);
                                else
                                    spy(D.dyst);
                                    pause(.01)
                                end
                            else
                                disp('d/dx matrix from s to t');
                                disp(full(D.dxst));
                                disp('d/dx matrix from t to s');
                                disp(full(D.dxts));
                                
                                disp('d/dy matrix from s to t');
                                disp(full(D.dyst));
                                disp('d/dy matrix from t to s');
                                disp(full(D.dyts));
                            end
%                             if BC1==BC2&&BC1==BC3&&BC1==BC4
%                                 disp(BC1);
%                                 pause(0.01);
%                             end
                        end
                    end
                
                elseif numDimensions ==1
                    disp([BC1 BC2]);
                    
                    myDomain = Domain(1,size1,[BC1, BC2],[0,0],2,1)
                    
                    disp('x vector:')
                    disp(myDomain.x)
                    D = DifferentialOperator(myDomain)
                    if visualization
                        figure;
                        spy(D.dxst);
                        pause(1);
                     
                    else
                        disp('d/dx matrix from s to t');
                        disp(full(D.dxst));
                        disp('d/dx matrix from t to s');
                        disp(full(D.dxts));
                    end
                end
            end
        end
    
end
