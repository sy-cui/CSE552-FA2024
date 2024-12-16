function A = init_acc(ModelA,NDOFT,neq,numnp,ndf)
    A = zeros(numnp, ndf);

    for node = 1:numnp
        for dir = 1:ndf
            gDOF = NDOFT(node, dir);
            if gDOF <= neq
                A(node,dir) = ModelA(gDOF);
            end
        end
    end
end