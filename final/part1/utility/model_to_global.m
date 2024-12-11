function u = model_to_global(model_x,model_dc,NDOFT,numnp,ndf,neq)
    % Map a variable in the equation indices to global indices
    u = zeros(numnp,ndf);
    for node = 1:numnp
        for dir = 1:ndf
            gDOF = NDOFT(node, dir);
            if gDOF <= neq
                u(node, dir) = model_x(gDOF,1);
            else
                u(node, dir) = model_dc(gDOF - neq);
            end
        end
    end
end