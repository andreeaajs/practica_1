def euler_explicit(dt, dx, f1, f0, f_1):
    return (dt/dx**2) * (f1 - 2*f0 + f_1) + dt 
