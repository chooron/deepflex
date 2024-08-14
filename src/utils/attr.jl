#* used for getting element attr
get_ode_func(::AbstractElement) = nothing
get_ode_func(ele::AbstractHydroElement) = ele.ode_func