function [ue,ElemFlag] = extract_elem_var(u, ie, ndm, nel, conn)
    % Extract element variables ('ue') from a global array ('u') based on
    % mesh connectivity ('conn')
    % 
    % Inputs:
    %
    %   u:          Global array containing all nodal information for this
    %               particular variable, e.g. NodeTable in FEA_Program.m.
    %
    %   ie:         Index of the element
    %   
    %   ndm:        Number of spatial dimensions
    %
    %   nel:        Number of nodes in this element
    %
    %   conn:       Mesh connectivity. See 'ix' in FEA_Program.m. 
    %
    % Outputs:
    %
    %   ue:         Subset of 'u' pertaining to the element of interest.
    %               Size (ndm x nel)
    %
    %   ElemFlag:   Index mapping between 'ue' and global array. For 
    %               example, k-th element of 'ElemFlag' contains the INDEX 
    %               corresponding to column of 'ue' in the global array 
    %               'u'. 
    %

    ue = zeros(ndm, nel);
    ElemFlag = zeros(nel, 1);

    for a = 1:nel
        node = conn(ie,a);
        ElemFlag(a) = node;
        for i = 1:ndm
            ue(i,a) = u(node,i);
        end
    end

end