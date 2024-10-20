# CEE 576 Homework 3
The source code for this homework implements the arc-length method for stepping along a scalar nonlinear function. 

## Helper functions

`fint.m`: Defines the nonlinear function $N(u)$. The function accepts vector inputs. Usage:
```m
u = linspace(0, 6, 128);
f = fint(u);
plot(u, f);
```

`dfint.m`: Defines the derivative (consistent tangent) of the function $N(u)$. The function accepts vector inputs. Usage:
```m
u = linspace(0, 6, 128);
df = dfint(u);
plot(u, df);
```

`arc.m`: Defines the arc-length function $f(\delta u, \delta \lambda)$. Requires inputs `Du` ($\delta u$), `Dl` ($\delta \lambda$), `dK` (diagonal of tangent matrix $k$), `b` and `c` (parameters in the arc-length function). Usage:
```m
% Suppose that da, K and q = K\F is defined 
dK = diag(K);
b = 0.5; 
c = (1 - b) / (q'*dK*q);
Dl = da; 
Du = da * q;
f = arc(Du,Dl,dK,b,c);
```

`darc.m`: Defines the gradient of arc-length function $f(\delta u, \delta \lambda)$ with respect to $\delta u$ and $\delta \lambda$. Requires inputs `Du` ($\delta u$), `Dl` ($\delta \lambda$), `dK` (diagonal of tangent matrix $k$), `b` and `c` (parameters in the arc-length function). Usage:
```m
% Suppose that da, K and q = K\F is defined 
dK = diag(K);
b = 0.5; 
c = (1 - b) / (q'*dK*q);
Dl = da; 
Du = da * q;
[dfDu, dfDl] = darc(Du,Dl,dK,b,c);
```

## Main solution source code
`hw3_p4.m`. The code implements the modified-Newton-Raphson method for the arc-length technique. 
Parameters are decorated with docstrings in file. 
Modifiable parameters are located at the top of file within the `Parameter` section. 
To run the script, in MATLAB console:
```sh
>>> hw3_p4.m
```
Output figure 1 contains the solution history, and figure 2 contains residual decay at each load step. 

