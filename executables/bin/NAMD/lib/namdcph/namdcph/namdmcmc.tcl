# namdmcmc.tcl
#
# Markov chain Monte Carlo and related statistics routines. 
#
# References:
# 1) Y. Chen and B. Roux "Efficient Hybrid Non-Equilibrium Molecular Dynamics -
#     Monte Carlo Simulations with Symmetric Momentum Reversal," J. Chem. Phys.
#     2014, 141, 114107.
#
source [file join [file dirname [info script]] "numtcl.tcl"]

# Check that this is being run through NAMD.
if {[catch numPes]} {
    puts "Expected a NAMD-Tcl interpreter - exiting."
    return 1
}
# =============================================================================
# Nonequilibrium switching routines
# =============================================================================
# runprpswitch
#
# Switch between states. With equal probability, either flip velocities before 
# _and_ after the switch or leave them unchanged.
#
# Arguments:
# ----------
# numsteps : int
#   Number of steps in the switch
#
# Returns:
# --------
# None 
# 
proc runprpswitch {numsteps} {
    set scale [expr {pow(-1, round(rand()))}]
    if {[expr {$scale < 0}]} {
        rescalevels $scale
        run $numsteps
        rescalevels $scale
    } else {
        run $numsteps
    }
}

# runfpswitch
#
# Switch between states and then flip the velocities.
#
# Arguments:
# ----------
# numsteps : int
#   Number of steps in the switch
#
# Returns:
# --------
# None 
# 
proc runfpswitch {numsteps} {
    run $numsteps
    rescalevels -1.
}

# =============================================================================
# Monte Carlo routines
# =============================================================================
# metropolisAcceptance
#
# Accept or reject a Monte Carlo move based on the change in energy using the 
# Metropolis criteria.
# 
# Arguments:
# ----------
# deltau : float
#   The _reduced_ change in energy (i.e. including inverse temperature 
#   factors).
#
# Returns:
# --------
# accept : bool
#   True/False value for accept/reject test 
# 
proc metropolisAcceptance {deltau} {
    if {$deltau <= 0.} {
        return 1
    } else {
        set Pacc [expr {exp(-$deltau)}]
        set Rand [expr {rand()}]
        return [expr {$Rand <= $Pacc}]
    }
}

# normalizeLogWeights
#
# Return a set of normalized probability weights given a set of _unnormalized_
# log weights. The stability of the calculation is aided by the introduction of
# an arbitrary shift factor, usually taken as the maximum of the set.
#
proc normalizeLogWeights {logWeights shift} {
    # Compute the (scaled) sum of the weights.
    set normedWeights [lrepeat [llength $logWeights] 0.0]
    set weightSum 0.0
    set i 0
    foreach logWeight $logWeights {
        set shiftedWeight [expr {exp($logWeight - $shift)}]
        lincr normedWeights $i $shiftedWeight
        set weightSum [expr {$weightSum + $shiftedWeight}]
        incr i
    }
    # Normalize the weights.
    lmultiply normedWeights [expr {1.0 / $weightSum}]
    return $normedWeights
}

