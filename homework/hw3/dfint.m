function dy = dfint(u)
    dy = 0.03*fint(u) + (0.57*u.^2-4.2*u+5.8)*exp(0.03*u);
end