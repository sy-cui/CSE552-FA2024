function [P,V,A] = predict_kinematics( ...
    pos,vel,acc,hht,NodeTable,model_dc,NDOFT,neq)

    P = pos;
    V = vel;
    A = acc;
    beta = hht(2); gamma = hht(3); dt = hht(4);
    
    numnp = size(NodeTable,1);
    ndf = size(NodeTable,2);

    for node = 1:numnp
        for dir = 1:ndf
            gDOF = NDOFT(node, dir);
            if gDOF <= neq
                P(node,dir) = pos(node,dir) + dt*vel(node,dir) + (0.5-beta)*dt^2*acc(node,dir);
                V(node,dir) = vel(node,dir) + dt*(1-gamma)*acc(node,dir);
                A(node,dir) = 0;
            else
                P(node,dir) = NodeTable(node,dir) + model_dc(gDOF - neq);
                V(node,dir) = 0;
                A(node,dir) = 0;
            end
        end
    end
end