function [fint_elem] = fint_elem_nh_2d(xc,xr,mateprop,nel,lint,framework)
    fint_elem = 0;

    for l = 1:lint
        if nel == 3
            [Wgt,r,s] =  intpntt(l,lint,0);
        else
            [Wgt,r,s] =  intpntq(l,lint,0);
        end

        if strcmpi(framework, 'tl')
            % Evaluate first derivatives of basis functions w.r.t reference
            [Qxy, Jfe] = shpg_2d(shpl_2d(r,s,nel),xr,nel);
            
            % Evaluate deformation gradient
            [F, J] = def_grad_2d(r,s,xc,xr);
    
            % Evaluate PK2 stress tensor
            S = pk2(F,J,mateprop);
    
            % Compute integrand at integration point
            f = reshape(F * S * Qxy(1:2,:), [], 1);

        elseif strcmpi(framework, 'ul')
            % Evaluate first derivatives of basis functions w.r.t current
            [Qxy, Jfe] = shpg_2d(shpl_2d(r,s,nel),xc,nel);

            % Evaluate deformation gradient
            [F, J] = def_grad_2d(r,s,xc,xr);
            
            % Compute Cauchy stress
            sig = cauchy(F,J,mateprop);

            % Compute integrand at integration point
            f = reshape(sig * Qxy(1:2,:),[],1);

        else
            error("'framework' must be one of 'ul' and 'tl'")
        end   

        fint_elem = fint_elem + Wgt*Jfe*mateprop(3)*f;

    end % quad point loop 

end
