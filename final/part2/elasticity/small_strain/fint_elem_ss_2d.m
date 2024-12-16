function fint_elem = fint_elem_ss_2d(xc,xr,mateprop,nel,lint)
    fint_elem = 0;

    u = reshape(xc - xr, [], 1);

    a = mateprop(1);
    b = mateprop(2);
    c = mateprop(3);

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

        strain = Bmat * u;
        eps11 = strain(1); 
        eps22 = strain(2); 
        eps12 = strain(3)/2;
        stress11 = a*(eps11+eps22)+b*eps11+c*(eps11^2+eps12^2);
        stress22 = a*(eps11+eps22)+b*eps22+c*(eps12^2+eps22^2);
        stress12 = b*eps12+c*eps12*(eps11+eps22);

        f = Bmat'*[stress11; stress22; stress12];
        
        fint_elem = fint_elem + Wgt*Jfe*t*f;

    end % quad point loop 

end