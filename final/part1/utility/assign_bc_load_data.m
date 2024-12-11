function [Fd,ModelDc,neq,nieq,NDOFT] = assign_bc_load_data( ...
    numnp, ndf, NodeBC, NodeLoad, BCLIndex)
    % Copyright (C) Arif Masud and Tim Truster
    %
    % 03/28/2009
    % 
    % Parse the boundary conditions (displacement and traction) and map
    % them to proper equation numberings. 
    % 
    % The final system of equation
    % should be of the form 
    %       | K_dd K_df | | d_1 | = | f_1 |
    %       | K_fd K_ff | | d_2 |   | f_2 |
    % where d_2 and f_1 can be determined from the Dirichlet and Neumann 
    % boundary conditions, respectively. The task of this routine is to 
    % find d_2 and f_1, presently named 'ModelDc' and 'Fd', based on the
    % boundary equations and rearrange them as the form above. 
    % 
    %
    % Inputs:
    % 
    %   numnp:      Total number of nodes in the mesh.
    %
    %   ndf:        Degrees of freedom in each node.
    % 
    %   NodeBC:     Dirichlet (displacement) boundary conditions formatted 
    %               as a matrix. Each row represents a Dirichlet condition,
    %               with three entires in the sequence [node dir value].
    %
    %   NodeLoad:   Neumann (traction) boudnary conditions formatted as a
    %               matrix. Each row represents a Neumann condition, with 
    %               three entires in the sequence [node dir value].
    %
    %   BCIndex:    Two element array. First entry is the number of
    %               Dirchlet boundaries (number of rows of 'NodeBC'), and
    %               the second entry is the number of Neumann boundary
    %               (number of rows in 'NodeLoad').
    % 
    % Outputs: 
    %   
    %   Fd:         RHS force vector of size (neq x 1) excluding the node-
    %               direction pairs with prescribed displacements.
    %               
    %   
    %   ModelDc:    Known displacements from the Dirichelt boundaries. 
    % 
    %   neq:        Number of unsolved equations.
    %
    %   nieq:       Number of "solved" equations for which the
    %               displacements are known but reaction forces are
    %               unknown.
    %   NDOFT:      A lookup-table that contains node-to-equation mapping.
    %               Each row represents a node consistent with global
    %               numbering, and the row contains 2*ndf entries. The
    %               first ndf entries represents mapping of node-direction
    %               pairs to equation indices, and last ndf entries
    %               contain the prescribed displacements, if applicable. 
    %               Example:
    %                   NDOFT = [6 7 1 | 2 0 0
    %                            ...   |   ...]
    %               Interpretation:
    %                   Node 1 direction 1 is mapped to equation 6.
    %                   Displacement is known to be 2. 
    %                   Node 1 direction 2 is mapped to equation 7. 
    %                   Displacement is known to be 0.
    %                   Node 1 direction 3 is mapped to equation 1. 
    %                   Displacement is unknown.
    %                   

    NDOFT = zeros(numnp,2*ndf);
    nieq = 0;
    
    % Get BC
    len = BCLIndex(1);
    for i = 1:len
        node = NodeBC(i,1);
        dir = NodeBC(i,2);
        displacement = NodeBC(i,3);
        NDOFT(node,dir) = -1;
        NDOFT(node,dir+ndf) = displacement;
        nieq = nieq + 1;
    end
    
    % Assign DOF
    neq = ndf*numnp - nieq;
    ModelDc = zeros(nieq, 1);
    na = 0;
    ni = neq;
    
    % Equations are numbered direction first: [x1 y1 x2 y2 ...]
    for i = 1:numnp
        for j = 1:ndf
            if NDOFT(i,j) == 0
                na = na + 1;
                NDOFT(i,j) = na;
            elseif NDOFT(i,j) == -1
                ni = ni + 1;
                NDOFT(i,j) = ni;
                ModelDc(ni-neq) = NDOFT(i,j+ndf);
            else
                % The current code should not reach here
                NDOFT(i,j) = NDOFT(NDOFT(i,j),j);
            end
        end
    end
    
    % Get Loads
    Fd = zeros(neq, 1);
    
    len = BCLIndex(2);
    for i = 1:len
        node = NodeLoad(i,1);
        dir = NodeLoad(i,2);
        force = NodeLoad(i,3);
        dof = NDOFT(node,dir);
        if dof <= neq
            Fd(dof) = force;
        end
    end

end