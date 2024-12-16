function u = update_mesh(NodeCurr,NodeTable,model_x,model_dc,NDOFT,neq)
    % Map a variable in the equation indices to global indices
    u = NodeCurr;
    numnp = size(NodeCurr,1);
    ndf = size(NodeCurr,2);
    for node = 1:numnp
        for dir = 1:ndf
            gDOF = NDOFT(node, dir);
            if gDOF <= neq
                u(node,dir) = NodeCurr(node,dir) + model_x(gDOF,1);
            else
                u(node,dir) = NodeTable(node,dir) + model_dc(gDOF - neq);
            end
        end
    end
end