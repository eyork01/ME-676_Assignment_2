# project2.jl
# Author: Ethan York
# University of Kentucky
# ME 676, Dr. Poonawala
# March 27, 2023

# Included Packages
using LinearAlgebra
using RigidBodyDynamics

# Global Variables
global q1 = [0.01,-0.5,-0.0,-2.0,-0.3,1.5,-0.7,0.1,0.1]
global q2 = [0.0,0.0,0.0,0.0,0.0,pi,0.01,0.01,0.01]
global zero_state = [0,0.0,0.0,0.0,0.0,0,0,0,0]

# Trajectory function
# Takes in a time and returns the desired position at that time
function traj(time)
    # Ensure that time is between 0 and 10
    time = max(min(time, 10), 0)
    # Current position = initial position + (% of Δq)
    return q1 .+ (time/10).*(q2-q1)
end

# PD Controller Function
function control_PD!(τ, t, state)
    # Compute a value for τ
    kp = 320
    kd = 20
    τ .= -kd .* velocity(state) - kp*(configuration(state) - traj(t))

    # Saturate
    act_sat = 50; # Actuator limits
    τ .= map( x -> x > act_sat ? act_sat : x,τ)
    τ .= map( x -> x < -act_sat ? -act_sat : x,τ)
end

# CTC Controller Function
function control_CTC!(τ, t, state)
    # Compute a value for τ
    kp = 550
    kd = 20
    PD = -kd .* velocity(state) - kp*(configuration(state) - traj(t))
    a = similar(velocity(state))
    a .= PD
    τ .= inverse_dynamics(state,a)

    # Saturate
    act_sat = 50; # Actuator limits
    τ .= map( x -> x > act_sat ? act_sat : x,τ)
    τ .= map( x -> x < -act_sat ? -act_sat : x,τ)
end

# MAIN PROGRAM
# Load mechanism and visualize at q0 config
delete!(vis)
mvis, mechanism = display_urdf("panda.urdf",vis)
state = MechanismState(mechanism)
set_configuration!(state, q1)
set_configuration!(mvis, configuration(state))

# Solve trajectory with PD controller
problemPD = ODEProblem(Dynamics(mechanism,control_PD!), state, (0.00,10.00))
solPD = transpose(solve(problemPD, Vern7()))

# Reset Configuration
set_configuration!(state, q1)
set_configuration!(mvis, configuration(state))

# Solve trajectory with CTC
problemCTC = ODEProblem(Dynamics(mechanism,control_CTC!), state, (0.00,10.00))
solCTC = transpose(solve(problemCTC, Vern7()))

# Print Solution
println("PD Norm: ", norm(solPD[end,1:9]-q2))
println("CTC Norm: ", norm(solCTC[end,1:9]-q2))