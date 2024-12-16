function fint_elem = fint_elem_lin_2d(xc,xr,mateprop,nel,lint,PSPS)
    fint_elem = 0;

    u = reshape(xc - xr, [], 1);

    E = mateprop(1);
    nu= mateprop(2);
    t = mateprop(3);

    if PSPS == 's' % Plane Stress
        Dmat = E/(1-nu^2)*[1  nu 0;
                           nu 1  0;
                           0  0  (1-nu)/2];
    
    else % Plane Strain
        Dmat = E/(1+nu)/(1-2*nu)*[1-nu nu    0;
                                  nu   1-nu  0;
                                  0    0     (1-2*nu)/2];     
    end

    for l = 1:lint
        if nel == 3
            [Wgt,r,s] =  intpntt(l,lint,0);
        else
            [Wgt,r,s] =  intpntq(l,lint,0);
        end

        [Qxy, Jfe] = shpg_2d(shpl_2d(r,s,nel),xr,nel);
        
        if nel == 3
            Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        
                    0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3)
                    Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3)];
        else
            Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        Qxy(1,4) 0 
                    0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3) 0        Qxy(2,4)
                    Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3) Qxy(2,4) Qxy(1,4)];
        end

        f = Bmat'*Dmat*Bmat*u;
        
        fint_elem = fint_elem + Wgt*Jfe*t*f;

    end % quad point loop 
end