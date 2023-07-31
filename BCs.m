classdef BCs
    %BCS Enum of possible boundary conditions
    %   Combinations of symmetric/antisymmetric boundaries, and whether
    %   the boundary is placed at an S- or T-point
    %   S corresponds to circle, T to triangle
    
    %   In 1D (fields transverse to boundaries):
    
    %   In the case of a symmetric boundary on the s point:
    
    %   The field corresponding to the s grid will be symmetric, or have 
    %   slope == 0 at the boundary (von Neumann). 
    %   This essentially amounts to no change, since it is assumed in the 
    %   boundary (outside the computational domain), the field is symmetric
    %   to the field in the domain, with the boundary as the symmetry point
    
    %   In this same case, the field corresponding to the t grid will be
    %   antisymmetric, or will have value forced to 0 (Dirichlet)
    %   Because this field is not defined at the boundary, to simulate
    %   forcing to 0, we imagine an equal negative value on the other side
    %   of the boundary when doing a derivative which takes us from t to s.
    
    %   In the case of an antisymmetric boundary on the s point:
    
    %   The field corresponding to the s grid will be antisymmetric,
    %   or will have value forced to 0 (Dirichlet).
    %   This simply corresponds to setting the boundary s point to 0.
    
    %   In this same case, the field corresponding to the t grid will be
    %   symmetric. This amounts to no change, since it is assumed the field
    %   on the other imagined side of the boundary is symmetric.
    
    %   In the case of a boundary on the t point:
    %   Everything above is the same, except s and t are reversed.
    
   enumeration
      antiSymS, symS, antiSymT, symT %Only implement with respect to s. 
   end
end

