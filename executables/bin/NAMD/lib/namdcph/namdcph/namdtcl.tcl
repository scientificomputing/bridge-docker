# namdtcl.tcl
#
# This file contains routines that mimic idiomatic NAMD behavior or otherwise
# re-implement utilities from the NAMD source code in pure Tcl.
#
# Check that this is being run through NAMD.
if {[catch numPes]} {
    puts "Expected a NAMD-Tcl interpreter - exiting."
    return 1
}

# =============================================================================
# Constants (These usually come from src/common.h)
# =============================================================================
variable BOLTZMANN 0.001987191 ;# in kcal/mol-K
variable PRESSUREFACTOR 69500.0 ;# convert kcal/mol to bar*A^3
variable PI 3.141592653589793
variable LN10 [expr {log(10)}]
variable LN2 [expr {log(2)}]

variable dummyXSC "# NAMD extended system configuration output file\n#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z\n0 0 0 0 0 0 0 0 0 0 0 0 0\n"

# =============================================================================
# Standard energy callback
# =============================================================================
global EnergyLabels
global EnergyValues
global energyArray

# energyCallback
#
# Set the "EnergyLabels" and "EnergyValues" lists to contain the labels and 
# values for the most recent NAMD MD step, respectively. In practice this is 
# only meant to be used as an argument to the NAMD built-in "callback" 
# procedure.
#
# Arguments:
# ----------
# labels : list of strings
#   a list of energy term labels (e.g. {TS BOND ANGLE ... ELECT ...})
# values : list of floats
#   a list of energy term values corresponding to the above labels
#
# Returns:
# --------
# None
#
proc energyCallback {labels values} {
    set ::EnergyLabels $labels
    set ::EnergyValues $values
    return
}

# storeEnergies
#
# Use the "EnergyLabels" and "EnergyValues" lists to populate the "energyArray"
# array. This is mostly a convenience function to cleanly and clearly retrieve
# energy information from NAMD. While NAMD will automatically call 
# energyCallback, this proc should be called each time "energyArray" is
# to be queried.
#
proc storeEnergies {} {
    foreach label $::EnergyLabels value $::EnergyValues {
        set ::energyArray($label) $value
    }
    return
}

# =============================================================================
# File handling
# =============================================================================
# Backup a file in the usual NAMD manner by appending ".BAK"
proc namdFileBackup {filename} {
    if {[file exists $filename]} {
        file copy -force $filename "$filename.BAK"
    }
    return
}

# =============================================================================
# Input checking (mimic functionality in src/SimParameters.C)
# =============================================================================
# Verify that the input value can be interpreted as a number. 
proc checkIsNumeric {name value} {
    set IsNumeric 0
    if {![catch {expr {abs($value)}}]} {
        set IsNumeric 1
    }
    set value [string trimleft $value 0]
    if {![catch {expr {abs($value)}}]} {
        set IsNumeric 1
    }
    if {!$IsNumeric} {
        abort "$name must be a numeric argument! (got $value)"
    }
    return 1
}

# Verify that the input value is a positive (nonzero) number.
proc checkIsPositive {name value} {
    checkIsNumeric $name $value
    if {$value <= 0} {
        abort "$name must be greater than zero!"
    }
    return 1
}

# Verify that the input value is a non-negative (possibly zero) number.
proc checkIsNotNegative {name value} {
    checkIsNumeric $name $value
    if {$value < 0} {
        abort "$name must be greater than or equal to zero!"
    }
    return 1
}

# Verify that a variable argument list is a multiple of the given number.
#
# Example: >> checkArglistIsMultiple argList 2
#          returns true if argList has an even number of elements.
#
proc checkArglistIsMultiple {argList multi} {
    if {[expr {[llength $argList] % $multi}] != 0} {
        abort "argument list must have multiple of $multi arguments!"
    }
    return 1
}

# =============================================================================
# Advanced parameter introspection
#
#   The following procs aim to solve the problem of keywords with redundant
# meaning. For example, multiple thermostats are available in NAMD and each has
# its own temperature variable. However, one generally does not care _how_ the
# temperature is regulated, but merely that it _is_ regulated and has a clear
# thermodynamic value. This is solved by having auxillary global tcl variables
# based on an exhaustive query of the keywords. For convenience, an extra 
# variable is also defined to give an explicit name to the thermostat so that
# consistent capitalization and formatting is achieved.
#
# Example: What is the temperature?
# >> # some MD operations
# >> getThermostat ;# populated the global variables
# >> set theCurrentTemperature [$::thermostatTempCmd] ;# query the temperature
# >> $::thermostatCmd off ;# turn off the termostat
#
# =============================================================================
global thermostatName ""
global thermostatCmd ""
global thermostatTempCmd ""

global barostatName ""
global barostatCmd ""
global barostatPresCmd ""
global barostatTempCmd ""

# getThermostat
#
# This only requires that the thermostat itself has been invoked, not that a
# temperature has been set. By default, the much maligned Berendsen thermostat
# raises an error, since this is generally incompatible with canonical
# sampling.
#
# Return 1 if a thermostat is set, else return 0.
#
proc getThermostat {{forbidBerendsen true}} {
    global thermostatName ""
    global thermostatCmd ""
    global thermostatTempCmd ""

    if {[isset langevin] && [langevin]} {
        set thermostatName "Langevin"
        set thermostatCmd langevin
        set thermostatTempCmd langevinTemp
    } elseif {[isset loweAndersen] && [loweAndersen]} {
        set thermostatName "Lowe-Andersen"
        set thermostatCmd loweAndersen
        set thermostatTempCmd loweAndersenTemp
    } elseif {[expr {[reassignFreq] > 0}]} {
        set thermostatName "Andersen (massive collisions)"
        set thermostatCmd reassignFreq
        set thermostatTempCmd reassignTemp
    } elseif {[isset tCouple] && [tCouple]} {
        if {$forbidBerendsen} {
            abort "Berendsen thermostat does not work for canonical sampling."
        } else {
            set thermostatName "Berendsen"
            set thermostatCmd tCouple
            set thermostatTempCmd tCoupleTemp
        }
    } else {
        return 0 
    }
    return 1
}

# getBarostat
#
# This only requires that the barostat itself has been invoked, not that a
# pressure has been set. By default, the slightly less maligned Berendsen
# barostat raises an error, since this is generally incompatible with
# isobaric-isothermal sampling.
#
# Return 1 if a barostat is set, else return 0.
#
proc getBarostat {{forbidBerendsen true}} {
    global barostatName ""
    global barostatCmd ""
    global barostatPresCmd ""
    global barostatTempCmd ""

    if {[isset langevinPiston] && [langevinPiston]} {
        set barostatName "Langevin Piston"
        set barostatCmd langevinPiston
        set barostatPresCmd langevinPistonTarget
        set barostatTempCmd langevinPistonTemp
    } elseif {[isset berendsenPressure] && [berendsenPressure]} {
        if {$forbidBerendsen} {
            abort "Berendsen barostat does not work for NpT sampling."
        } else {
            set barostatName "Berendsen"
            set barostatCmd berendsenPressure
            set barostatPresCmd berendsenPressureTarget
        }
    } else {
        return 0 
    }
    return 1
}

# Alchemical Interactions
# 
#   NAMD uses separate scaling schemes for bonded, electrostatic, and van der 
# Waals terms, each of which is determined by a separate parameter. This
# essentially permits staggered scaling so that some interactions can be turned
# off quickly (i.e. electrostatics) while others can be turned off more slowly
# (i.e. bonded and van der Waals). The following routines just re-implement 
# the NAMD internals so that only a single scaling parameter (alchLambda) needs
# to be tracked.
#
#   For convenience, all 6 scaling factors can also be obtained at the same 
# time - this is a bit opaque sometimes, but probably much cleaner than using
# separate variables for each energy component. 
#   There are 3 interaction types  (bonds, electrostatics, van der Waals - in 
# _that_ order) and two alchemical groups (1 and 2 - in _that_ order).
#
proc getLambdas {lambda1} {
    set lambda2 [expr {1. - $lambda1}]
    return [list [getBondLambda $lambda1] [getElecLambda $lambda1]\
                 [getVdwLambda $lambda1] [getBondLambda $lambda2]\
                 [getElecLambda $lambda2] [getVdwLambda $lambda2]\
           ]
}

proc getElecLambda {lambda} {
    set ls [expr [alchElecLambdaStart]]
    return [expr {($lambda <= $ls) ? 0. : [expr {($lambda-$ls) / (1.-$ls)}]}]
}

proc getVdwLambda {lambda} {
    set le [expr [alchVdwLambdaEnd]]
    return [expr {($lambda >= $le) ? 1. : [expr {$lambda / $le}]}]
}

proc getBondLambda {lambda} {
    set le [expr [alchBondLambdaEnd]]
    return [expr {($lambda >= $le) ? 1. : [expr {$lambda / $le}]}]
}

# Return the alchemical energy components in the same list structure used by
# getLambdas. Note that energyArray should be passed _by name_, not by value.
#
# Example:
# callback energyCallback
# storeEnergies ;# this populates the global array energyArray
# set energyList [getAlchEnergies ::energyArray]
# 
proc getAlchEnergies {energyArray} {
    upvar 1 $energyArray energies
    return [list $energies(BOND1) $energies(ELEC1) $energies(VDW1)\
                 $energies(BOND2) $energies(ELEC2) $energies(VDW2)\
           ]
}

