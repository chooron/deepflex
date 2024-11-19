using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots

# Define the system with discrete updates
function discrete_condition(u, t, integrator)
    # Trigger condition for discrete updates
    return true
end

function discrete_affect!(integrator)
    # Get current values
    Q = integrator.u[1]
    I = integrator.p[4]*Q + integrator.p[5]  # A*Q + R
    
    # Calculate next Q using discrete equation
    C0, C1, C2, A, R = integrator.p
    Q_next = C0*(A*Q + R) + C1*I + C2*Q
    
    # Update state
    integrator.u[1] = Q_next
end

# Create the discrete callback
cb = DiscreteCallback(discrete_condition, discrete_affect!, save_positions=(false,true))

# Define the continuous system (simplified for discrete updates)
function continuous_system!(du, u, p, t)
    du[1] = 0.0  # No continuous evolution between discrete steps
end

# Parameters and initial conditions
p = (0.2, 0.3, 0.4, 1.0, 2.0)  # C0, C1, C2, A, R
u0 = [1.0]  # Initial Q
tspan = (0.0, 10.0)

# Create and solve the problem with discrete callbacks
prob = ODEProblem(continuous_system!, u0, tspan, p)
sol = solve(prob, Tsit5(), callback=cb, tstops=0:1:10)

# Calculate I values for plotting
I_values = [p[4]*q[1] + p[5] for q in sol.u]
t_values = sol.t

# Plot results
plot(t_values, [first.(sol.u) I_values], 
     label=["Q(t)" "I(t)"],
     title="Discrete System Solution",
     xlabel="Time",
     ylabel="Value",
     marker=:circle)