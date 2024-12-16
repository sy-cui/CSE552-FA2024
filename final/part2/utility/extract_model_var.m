function [u_eq, u_ieq] = extract_model_var(field,NDOFT,neq)
    numnp = size(field,1);
    ndf = size(field,2);
    nst = numnp*ndf;

    u_eq = zeros(neq,1);
    u_ieq = zeros(nst-neq,1);

    for node = 1:numnp
        for dir = 1:ndf
            gDOF = NDOFT(node, dir);
            if gDOF <= neq
                u_eq(gDOF) = field(node,dir);
            else
                u_ieq(gDOF-neq) = field(node,dir);
            end
        end
    end
end