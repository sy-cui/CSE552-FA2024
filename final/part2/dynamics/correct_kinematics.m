function [P,V,A] = correct_kinematics( ...
    pos,vel,acc,hht,NodeTable,model_da,model_dc,NDOFT,neq)

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
                da = model_da(gDOF,1);
                P(node,dir) = pos(node,dir) + beta*dt^2*da;
                V(node,dir) = vel(node,dir) + gamma*dt*da;
                A(node,dir) = acc(node,dir) + da;
            else
                P(node,dir) = NodeTable(node,dir) + model_dc(gDOF - neq);
                V(node,dir) = 0;
                A(node,dir) = 0;
            end
        end
    end
end