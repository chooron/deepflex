using DataInterpolations

u = rand(5)
t = 0:4
interp = LinearInterpolation(u, t)
interp(3.5) # Gives the linear interpolation value at t=3.5
interp[3] # Gives the linear interpolation value at t=3.5