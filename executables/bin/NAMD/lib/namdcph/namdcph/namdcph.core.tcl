# namdcph.core.tcl
# 
# Core utilities for setting up and running constant pH simulations in NAMD.
#
package require Tcl 8.5

source [file join [file dirname [info script]] "namdtcl.tcl"]
source [file join [file dirname [info script]] "namdmcmc.tcl"]
source [file join [file dirname [info script]] "cphtitrator.tcl"]
source [file join [file dirname [info script]] "json.tcl"]

namespace eval ::namdcph {
    namespace import ::cphTitrator::*

    variable TITLE "namdcph)"
    variable SystempH
    variable configFilename ""
    variable restartFilename ""
    variable restartFreq 0
    variable outFile ""
    variable numMinSteps 0
    variable excludeList [list]
    variable AlchFrcCons 100.0
    variable MDBasename namdcph.md
    variable SWBasename namdcph.sw
    # cphSystem data and defaults
    variable stateInfo [dict create]
    # cphTitrator data and defaults
    variable moveInfo [dict create]
    dict set moveInfo maxProposalAttempts 0
    dict set moveInfo default numsteps 20 
    dict set moveInfo default weight 1.0

    variable residueAliases [dict create]
    dict set residueAliases HIS {HSD HSE HSP}

    namespace export *
}

# =============================================================================
# Main NAMD Routines
# =============================================================================
# ::namdcph::cphRun
#
# Run a set number of neMD/MC cycles.
#
proc ::namdcph::cphRun {numsteps {numcycles 1}} {
    # Initialize NAMD and build a constant pH enabled PSF.
    initialize
    if {$::namdcph::numMinSteps > 0} {
        set tmp [firstTimeStep]
        minimize $::namdcph::numMinSteps
        firstTimeStep $tmp
    }
    set cphlog [openCpHLog]
    set firstCycle 1
    set lastCycle [expr {$firstCycle + $numcycles - 1}] 
    for {set cycle $firstCycle} {$cycle <= $lastCycle} {incr cycle} { 
        # (1) Propose a move from the full move set.
        #
        lassign [proposeMove] accept swNumsteps segresidList nattempts
        # (2) At this point we have either selected a switch or rejected a
        # whole bunch of moves. If the former, perform the switch.
        #
        if {!$accept} {
            cphPrint "All proposals rejected ($nattempts total)."
            set runArgs [list norepeat $numsteps]
        } else {
#            printProposalSummary $segresidList
            cphPrint "Proposal accepted ($nattempts), attemping a switch."
            set accept [runSwitch $swNumsteps $segresidList]
            set runArgs [list $numsteps]
            # Only accumulate statistics for attempted switches.
            accumulateAcceptanceRate $accept $segresidList
        }
        cphSystem update $accept $segresidList
        # (3) Perform whatever equilibrium sampling is desired.
        #
        runMD {*}$runArgs
        # (4) Output cycle information and a checkpoint file if necessary.
        #
        puts $cphlog "[format "%6d" $cycle] [cphSystem get occupancy]"
        flush $cphlog
        writeRestart "[outputname].cphrst" $cycle
    }
    writeRestart force "[outputname].cphrst" $cycle
    # Cleanup temporary files
    file delete {*}[glob [getSWBasename].*]
    file delete {*}[glob [getMDBasename].*]
    close $cphlog
    printTitratorReport
    return
}

#
# FOR ADVANCED USE ONLY!!!
#
# ::namdcph::testResidue
#
# Test one or more residue definitions.
#
proc ::namdcph::testResidue {resnames {verbose 0}} {
    pH 7.0
    cphForceConstant 0.0
    # Initialize NAMD and build a constant pH enabled PSF.
    initialize
    $::thermostatTempCmd 0.0
    outputEnergies 1
    alchLambdaFreq 0
    set report [dict create]
    set terms [list BOND ANGLE DIHED IMPRP ELECT VDW POTENTIAL]
    foreach resLabel [cphSystem get reslabel] {
        lassign [split $resLabel ":"] segid resid residue
        if {[lsearch $resnames $residue] < 0} {
            continue
        } 
        set segresid "$segid:$resid"
        set states [cphSystem get stateList $segresid]
        foreach state0 $states {
            foreach state1 $states {
                if {[string match $state0 $state1]} {
                    continue
                }
                # Force the initial state as state0
                output [getMDBasename]
                psfgenRead [getMDBasename]
                cphSystem initializeState $segresid $state0
                cphSystem alchemifypsf $segresid 0.0 0.0
                cphSystem dealchemifypsf $segresid
                regenerate angles dihedrals 
                cphSystem update 1 $segresid
                psfgenWrite [getMDBasename] 
                reloadAndReinit [getMDBasename] false
                # Now perform an alchemifying and dealchemifying transition.
                run 0
                storeEnergies
                array set MDEnergies0 [array get ::energyArray]
                cphSystem set trialState $segresid $state1
                alchemify $segresid
                reloadAndReinit [getSWBasename] true
                alchLambda 0.0
                run 0
                storeEnergies
                array set SWEnergies0 [array get ::energyArray]
                alchLambda 1.0
                run 0
                storeEnergies
                array set SWEnergies1 [array get ::energyArray]
                dealchemify $segresid
                cphSystem update 1 $segresid
                run 0
                storeEnergies
                array set MDEnergies1 [array get ::energyArray]
                # Record the differences in energy.
                set label [format "%s%sH%s" $residue $state0 $state1]
                dict set report $label [dict create]
                dict set report $label 0 [dict create]
                dict set report $label 1 [dict create]
                foreach term $terms {
                    set diff0 [expr {$MDEnergies0($term)-$SWEnergies0($term)}]
                    set diff1 [expr {$MDEnergies1($term)-$SWEnergies1($term)}]
                    if {[expr {abs($diff0)}] < 1e-5} {
                        set diff0 0.0
                    }
                    if {[expr {abs($diff1)}] < 1e-5} {
                        set diff0 0.0
                    }
                    dict set report $label 0 $term $diff0
                    dict set report $label 1 $term $diff1
                }
            }
        }
    }
    dict for {label data} $report {
        cphPrint [format "%-14s  %14s" $label "RMS-error"]
        foreach lambda {0 1} label {alchemify dealchemify} {
            set i 1
            set termStr " "
            set diffStr " "
            set rmse 0.0
            dict for {term diff} [dict get $data 0] {
                set termStr [format "%s % 14s" $termStr $term]
                set diffStr [format "%s % 14.4f" $diffStr $diff]
                set rmse [expr {$rmse + $diff*$diff}]
                incr i
                if {$i == 5} {
                    set termStr "$termStr     "
                    set diffStr "$diffStr     "
                    set i 0
                }
            }
            set rmse [expr {2*$rmse / [dict size $data]}]
            cphPrint [format "%-14s: % 14.4f" $label $rmse]
            if {$verbose} {
                cphPrint $termStr
                cphPrint $diffStr
            }
        }
    }
    # Cleanup temporary files
    file delete {*}[glob [getSWBasename].*]
    file delete {*}[glob [getMDBasename].*]
    return 
}

# =============================================================================
# "Setter" Routines - used as new keywords in NAMD
#
#   All of these procedures take the desired value as an argument and return 
# that value. For default values see the decalaration in the namespace 
# header.
# =============================================================================
# -----------------
# Required Keywords
# -----------------
# ::namdcph::pH
#
# pH value for the simulation
#
proc ::namdcph::pH {pHValue} {
    checkIsPositive "pH" $pHValue
    variable ::namdcph::SystempH $pHValue 
    return
}

# ::namdcph::cphConfigFile
#
# Configuration file for constant pH residues.
#
proc ::namdcph::cphConfigFile {filename} {
    variable ::namdcph::configFilename [string trim $filename]
    return
}

# ::namdcph::cphNumstepsPerSwitch
#
# For odd arguments, the first argument is assumed to be a default switch time.
# All remaining arguments are presumed to be label/numsteps pairs for specific
# moves.
#
proc ::namdcph::cphNumstepsPerSwitch {args} {
    variable ::namdcph::moveInfo
    if {[expr {[llength $args] % 2}]} {
        set numsteps [lindex $args 0]
        checkIsNotNegative numSwitchSteps $numsteps
        dict set moveInfo default numsteps [expr {int($numsteps)}]
        set args [lrange $args 1 end]
    }
    checkArglistIsMultiple $args 2
    for {set i 0} {$i < [llength $args]} {incr i 2} {
        lassign [lrange $args $i [expr {$i+1}]] moveLabel numsteps
        checkIsNotNegative numSwitchSteps $numsteps
        dict set moveInfo $moveLabel numsteps [expr {int($numsteps)}]
    }
    return
}

# ---------------------
# Commonly Used Options
# ---------------------
# ::namdcph::cphSetResidueState
#
# Set the state of one or more residues using psfgen syntax.
#
proc ::namdcph::cphSetResidueState {args} {
    checkArglistIsMultiple $args 2
    variable ::namdcph::stateInfo
    for {set i 0} {$i < [llength $args]} {incr i 2} {
        lassign [lrange $args $i [expr {$i+1}]] segresid state
        dict set stateInfo $segresid state $state
    }
    return
}

# ::namdcph::cphSetResiduepKai
#
# Set the inherent pKa of one or more residues using psfgen syntax.
#
proc ::namdcph::cphSetResiduepKai {args} {
    checkArglistIsMultiple $args 2
    variable ::namdcph::stateInfo
    for {set i 0} {$i < [llength $args]} {incr i 2} {
        #NB pKai may be a list of values - check individually.
        lassign [lrange $args $i [expr {$i+1}]] segresid pKai
        foreach pKa $pKai {
            checkIsPositive pKai $pKa
        }
        dict set stateInfo $segresid pKai $pKai
    }
    return
}

# ::namdcph::cphRestartFile 
#
# Restart a constant pH run from a restart file.
#
proc ::namdcph::cphRestartFile {filename} {
    variable ::namdcph::restartFilename $filename
    return
}

# ::namdcph::cphRestartFreq
#
# Frequency (in neMD/MC cycles) at which to save constant pH restart files
#
proc ::namdcph::cphRestartFreq {frequency} {
    checkIsNotNegative cphRestartFreq $frequency
    variable ::namdcph::restartFreq $frequency
    return
}

# ::namdcph::cphOutFile
#
# Name for constant pH output file - default is [outputname].cphlog
#
proc ::namdcph::cphOutFile {filename} {
    variable ::namdcph::outFile $filename
    return
}

# ::namdcph::cphProposalWeight
#
# The (unnormalized) proposal weight assigned to each move. The default is for
# all such weights to be equal (uniformally distributed moves).
#
proc ::namdcph::cphProposalWeight {args} {
    variable ::namdcph::moveInfo
    checkArglistIsMultiple $args 2
    for {set i 0} {$i < [llength $args]} {incr i 2} {
        lassign [lrange $args $i [expr {$i+1}]] moveLabel weight
        checkIsNotNegative proposalWeight $weight
        dict set moveInfo $moveLabel weight [expr {1.*$weight}]
    }
    return
}

# ::namdcph::cphMaxProposalAttempts
#
# Number of attempted MC proposals from each move set before giving up.
#
# Values less than 1 default to the number of residues in the system.
# NB: This is not the necessarily the same thing as attempting a move once for
#     each residue.
#
proc ::namdcph::cphMaxProposalAttempts {maxAttempts} {
    checkIsNumeric cphMaxProposalAttempts $maxAttempts
    variable ::namdcph::moveInfo
    dict set moveInfo maxProposalAttempts [expr {int($maxAttempts)}]
    return
}

# ::namdcph::cphNumMinSteps
#
# Number of minimization steps to perform before dynamics.
#
# This is especially useful when the initial states are randomized according to
# the pH.
#
proc ::namdcph::cphNumMinSteps {numsteps} {
    checkIsNotNegative cphNumMinSteps $numsteps
    variable ::namdcph::numMinSteps $numsteps
    return
}

proc ::namdcph::cphExcludeResidue {args} {
    variable ::namdcph::excludeList
    foreach segresid $args {
         lappend excludeList $segresid
    }
    return
}

# -------------------
# Specialized Options
# -------------------
# ::namdcph::cphForceConstant
#
# Force constant for zero-length bonds between alchemical atoms in kcal/mol-A^2
#
proc ::namdcph::cphForceConstant {forceConstant} {
    checkIsNotNegative cphForceConstant $forceConstant
    variable ::namdcph::AlchFrcCons $forceConstant
    return
}

# ::namdcph::cphMDBasename
#
# Basename for (temporary) constant pH MD input files - the only output of
# interest to the user should be in the usual places
#
proc ::namdcph::cphMDBasename {basename} {
    variable ::namdcph::MDBasename $basename
    return
}

# ::namdcph::cphSwitchBasename
#
# Basename for (temporary) constant pH switch trajectory output files
#
proc ::namdcph::cphSwitchBasename {basename} {
    variable ::namdcph::SWBasename $basename
    return
}

# =============================================================================
# Nonequilibrium Switching Routines
# =============================================================================
# ::namdcph::runMD
proc ::namdcph::runMD {args} {
   run {*}$args 
}

# ::namdcph::runSwitch
#
# Run a neMD/MC switch on the specified residue(s).
#
# This includes multiple steps:
# (1) Modify I/O settings to differentiate from regular MD
# (2) Build alchemically enabled side chains for the desired residues
# (3) Perform a switching trajectory with appropriate momentum reversals
# (4) Perform Metropolis MC using an appropriate work quantity
# (5) Reinitialize MD based on the MC result
#
# Arguments:
# ----------
# numsteps : int
#   Number of steps in the switch
# segresidList : list 
#   One or more "<segid>:<resid>" specifications - this is the same syntax as 
#   for the regular psfgen patch command.
#
# Returns:
# --------
# accept : boolean
#   Result of MC accept/reject test 
#
proc ::namdcph::runSwitch {numsteps segresidList} {
    # (1) Checkpoint and modify output parameters. 
    #
    checkpoint
    storeEnergies
    set storedOutputEnergies [outputEnergies]
    set storedDCDFreq [dcdFreq]
    set storedTimestep $::energyArray(TS)
    outputEnergies $numsteps
    dcdFreq 0
    firstTimestep 0
    # (2) Build the alchemical switch inputs.
    #
    alchemify $segresidList
    # (3) Run the switch trajectory with momentum reversal.
    #
    runprpswitch $numsteps
    outputEnergies $storedOutputEnergies
    dcdFreq $storedDCDFreq
    # (4) Compute the work with state dependent energy shifts.
    #
    storeEnergies
    set DeltaE [cphSystem compute switch $segresidList]
    set Work [expr {$::energyArray(CUMALCHWORK) + $DeltaE}]
    set ReducedWork [expr {$Work / ($::BOLTZMANN*[$::thermostatTempCmd])}]
    set accept [metropolisAcceptance $ReducedWork] 
#    printProposalSummary $segresidList
    set tmp [printProposalSummary $segresidList]
    cphPrint [format "%s WorkCorr % 10.4f CorrWork % 10.4f"\
            [join $tmp "/"] $DeltaE $Work]
    # (5) Reinitialize for the next MD step based on accept/reject.
    #
    outputEnergies $storedOutputEnergies
    dcdFreq $storedDCDFreq
    firstTimestep $storedTimestep
    if {$accept} {
        cphPrint "Switch accepted!"
        dealchemify $segresidList
    } else {
        cphPrint "Switch rejected!"
        alch off
        revert
        reloadAndReinit [getMDBasename] false
    }
    return $accept
}

# ::namdcph::alchemify
#
# Reinitialize the system with alchemical sidechains for the given residues.
#
# This includes multiple steps:
# (1) The NAMD state is written to disk and re-read by PSFGEN.
# (2) Appropriate alchemical patches are applied and alchemical atom 
#     coordinates and velocities are sampled.
# (3) New NAMD inputs are written and then re-read. This includes new 
#     extraBonds for the alchemical sidechain.
#
# Arguments:
# ----------
# segresidList : list of strings
#   One or more "<segid>:<resid>" specifications - this is the same syntax as 
#   for the regular psfgen patch command.
#
# Returns:
# --------
# None
#
proc ::namdcph::alchemify {segresidList} {  
    set T [$::thermostatTempCmd]
    set FrcCons [getAlchFrcCons]
    alch on
    # (1) Read in the nonalchemical PSF and apply patches.
    #
    output [getMDBasename]
    psfgenRead [getMDBasename]
    # (2) Apply patches and build coordinates and velocities.
    #
    foreach segresid $segresidList {
        cphSystem alchemifypsf $segresid $FrcCons $T
    }
    regenerate angles dihedrals
    # (3) Write a new set of inputs and reinitialize.
    #
    psfgenWrite [getSWBasename] [mdFilename xsc]
    # Dummy atoms have been built and a PDB has been written. We can now query
    # atom indices and build extraBonds. If psfgenWrite were _not_ called, then
    # the atomid queries would all return zero (that would be bad).
    #
    set ExtraBondsFile [open [swFilename extrabonds] "w"]
    foreach segresid $segresidList {
        puts $ExtraBondsFile [cphSystem get alchBonds $segresid $FrcCons]
    }
    close $ExtraBondsFile
    reloadAndReinit [getSWBasename] true
    return
}

# ::namdcph::dealchemify
#
# Remove alchemical side
#
# Arguments:
# ----------
# segresidList : list of strings
#   One or more "<segid>:<resid>" specifications - this is the same syntax as 
#   for the regular psfgen patch command.
#
# Returns:
# --------
# None
#
proc ::namdcph::dealchemify {segresidList} {
    output [getSWBasename]
    psfgenRead [getSWBasename]
    foreach segresid $segresidList {
        cphSystem dealchemifypsf $segresid
    }
    psfgenWrite [getMDBasename] [swFilename xsc]
    alch off
    reloadAndReinit [getMDBasename] false
    return
}

# =============================================================================
# Constant pH specific I/O
# =============================================================================
# ::namdcph::readRestart
#
# Read in a constant pH state from file.
#
# Arguments:
# ----------
# restartFilename : string
#   Name of (JSON format) restart file to read
#
# Returns:
# --------
# stateList : list
#   State codes for each titratable residue
# pKaiList : list
#   Latest estimate of the inherent pKa
#
proc ::namdcph::readRestart {restartFilename} {
    set RestartFile [open $restartFilename "r"]
    set Restart [json::json2dict [read $RestartFile]]
    close $RestartFile
    if {[dict exists $Restart states]} {
        set stateList [dict get $Restart states]
    } else {
        set stateList {}
    }
    if {[dict exists $Restart pKais]} {
        set pKaiList [dict get $Restart pKais]
    } else {
        set pKaiList {}
    }
    if {[llength $stateList] > 0 && [llength $pKaiList] > 0} {
        if {[llength $stateList] != [llength $pKaiList]} {
            cphAbort "mismatch in states/pKais in $restartFilename"
        }
    }
    if {[dict exists $Restart MCmoves]} {
        set MCmoves [dict get $Restart MCmoves]    
    } else {
        set MCmoves [dict create]
    }
    return [list $stateList $pKaiList $MCmoves]
}

# ::namdcph::writeRestart
#
# Write the current constant pH state information to restartFilename based on 
# the restart frequency. If the "force" keyword precedes the arguments, ignore 
# the restart frequency and write no matter what.
#
# Arguments:
# ----------
# restartFilename : string
#   Name of (JSON format) restart file to write
# cycle : int
#   Last neMD/MC cycle index
#
# Returns:
# --------
# None
#
proc ::namdcph::writeRestart {args} { 
    variable ::namdcph::restartFreq
    if {[string match [lindex $args 0] force]} {
        set restartFilename [lindex $args 1]
        set cycle [expr {int([lindex $args 2])}]
    } else {
        set restartFilename [lindex $args 0]
        set cycle [expr {int([lindex $args 1])}]
        if {!($restartFreq) || ![expr {$cycle % $restartFreq}]} {
            return
        }
    }
    namdFileBackup $restartFilename
    set RestartFile [open $restartFilename "w"]
    # TODO: This is god-awful. Is there not a useable json encoder that can do
    # this for us with some guarantee of accuracy?
    set cycleStr "\"cycle\":$cycle"
    set stateStr "\"states\":\["
    foreach state [cphSystem get state] {
        set stateStr "$stateStr\"$state\","
    }
    set stateStr "[string trimright $stateStr ","]\]"
    set pKaiStr "\"pKais\":\["
    foreach pKaiList [cphSystem get pKai] {
        set pKaiStr "$pKaiStr\["
        foreach pKai $pKaiList {
            set pKaiStr "$pKaiStr$pKai,"
        }
        set pKaiStr "[string trimright $pKaiStr ","]\],"
    }
    set pKaiStr "[string trimright $pKaiStr ","]\]"
    set MCStr "\"MCmoves\":\{"
    dict for {moveLabel data} [getMoveSet] {
        set numsteps [dict get $data numsteps]
        set weight [dict get $data weight]
        set MCStr "$MCStr\"$moveLabel\":\{\"numsteps\":$numsteps,\"weight\":$weight\},"
    }
    set MCStr "[string trimright $MCStr ","]\}"

    puts $RestartFile "\{$cycleStr,$stateStr,$pKaiStr,$MCStr\}"
    close $RestartFile
    # Always checkpoint the PSF and PDB when a restart is written.
    file copy -force [mdFilename psf] "[outputName].psf"
    file copy -force [mdFilename pdb] "[outputName].pdb"
    return
}

# ::namdcph::openCpHLog
#
# Open a new constant pH log for proton occupancies and return the file object.
#
proc ::namdcph::openCpHLog {} {
    variable ::namdcph::outFile
    if {[string length $outFile] > 0} {
        set logFilename $outFile
    } else {
        set logFilename "[outputname].cphlog"
    }
    namdFileBackup $logFilename
    set cphlog [open $logFilename "w"]
    puts $cphlog "#pH [getpH]"
    puts $cphlog "#[join [cphSystem get reslabel] " "]"
    return $cphlog 
}

# =============================================================================
# Convenience Routines 
# =============================================================================
# ::namdcph::swFilename
#
# Get an appropriate filename for a temporary file used during switches
#
proc ::namdcph::swFilename {ext} {
    return "[getSWBasename].$ext"
}

# ::namdcph::mdFilename
#
# Get an appropriate filename for a temporary file used to launch regular MD
#
proc ::namdcph::mdFilename {ext} {
    return "[getMDBasename].$ext"
}

# ::namdcph::clearExtraBonds
#  
# Clear all bonds in the extraBondsFile.
#
proc ::namdcph::clearExtraBonds {} {
    set ExtraBondsFilename [swFilename extrabonds]
    cphPrint "clearing extraBondsFile $ExtraBondsFilename"
    set ExtraBondsFile [open $ExtraBondsFilename "w"]
    puts $ExtraBondsFile ""
    close $ExtraBondsFile
    return
}

# ::namdcph::reloadAndReinit
#
# Read in a new PSF/PDB pair and then reinitialize atoms using a basename read.
#
proc ::namdcph::reloadAndReinit {basename keepExtraBonds} {
    if {!$keepExtraBonds} {
        clearExtraBonds
    }
    structure "$basename.psf" pdb "$basename.pdb"
    reinitatoms $basename
    return
}

proc ::namdcph::cphPrint {text} {
    print "$::namdcph::TITLE $text"
}

proc ::namdcph::cphAbort {msg} {
    abort "$::namdcph::TITLE $msg"
}

# =============================================================================
# Setup Routines
# =============================================================================
# ::namdcph::checkSettings
#
# Check the settings from the configuration file and verify that they are
# compatible with constant pH.
#
proc ::namdcph::checkSettings {} {
    variable ::namdcph::configFilename
    # Check constant pH specific settings. 
    if {[string length configFilename] <= 0} {
        cphAbort "A constant pH configuration file is required."
    }
    if {![info exists ::namdcph::SystempH]} {
        cphAbort "A pH value is required."
    }

    getThermostat
    if {[getBarostat]} { ; # For now, disable constant pressure.
        cphAbort "Constant pH does not currently support the use of a barostat."
    }
    return
}

# ::namdcph::initialize
#
# Initialize the system for constant pH. This requires two main things to
# happen:
#
# 1) nonequilibrium alchemical transformations must be enabled
# 2) the PSF/PDB must be rebuilt to include dummy atoms and possibly reassigned
#    protonation states
#
proc ::namdcph::initialize {} {
    variable ::namdcph::configFilename
    variable ::namdcph::restartFilename
    variable ::namdcph::stateInfo
    variable ::namdcph::moveInfo
    variable ::namdcph::residueAliases
    variable ::namdcph::excludeList
    checkSettings
    callback energyCallback    
    # 1) Set up alchemical keywords and run an energy evaluation.
    #
    initializeAlch
    # 2) Rebuild the PSF with dummy protons and modify protonation states as 
    # needed. Build the residue definitions, assign states to each residue, and 
    # rebuild the topology to reflect those states.
    #
    cphPrint "initializing constant pH PSF..."
    if {[string length $restartFilename] > 0} {
        lassign [readRestart $restartFilename] states pKais MCmoves
        lassign [list 0.0 false] temp buildH
    } else {
        lassign [list {} {} [dict create]] states pKais MCmoves
        lassign [list [$::thermostatTempCmd] true] temp buildH
    }
    set ConfigFile [open $configFilename "r"]
    set templateDefs [json::json2dict [read $ConfigFile]]
    close $ConfigFile
    cphSystem build $templateDefs $residueAliases $excludeList
    # Use the system topology to make assignments from the restart
    if {[llength $states] && [llength $pKais]} { 
        foreach segresid [cphSystem get segresids] state $states pKai $pKais {
            if {![dict exists $stateInfo $segresid]} {
                dict set stateInfo $segresid state $state
                dict set stateinfo $segresid pKai $pKai
            } else {
                if {![dict exists $stateInfo $segresid state]} {
                    dict set stateInfo $segresid state $state
                }
                if {![dict exists $stateInfo $segresid pKai]} {
                    dict set stateInfo $segresid pKai $pKai
                }
            }
        }
    }
    cphSystem initialize [getpH] $temp $buildH $stateInfo
    # 3) Build the MC move set (the "titrator").
    set moveInfo [dict merge $MCmoves $moveInfo]
    set unusedInfo [buildTitrator [getpH] $moveInfo]
    # 4) Report to stdout.
    printSettingsSummary
    printSystemSummary
    printTitratorSummary
    if {[dict size [dict remove $unusedInfo default]]} {
        cphPrint "WARNING! Unused/unrecognized MC move info!"
        cphPrint [format "%-14s : %-s" "Move Label" "Info"]
        dict for {moveLabel data} [dict remove $unusedInfo default] {
            cphPrint [format "%-14s : %-s" $moveLabel $data]
        }
    }
    # 5) Write to disk and prepare for MD.
    if {[isset extendedSystem]} {
        set inputXSCName [extendedSystem]
    } else {
        set inputXSCName [mdFilename xsc]
        set inputXSC [open $inputXSCName "w"]
        puts $inputXSC $::dummyXSC 
        close $inputXSC
    }
    psfgenWrite [getMDBasename] $inputXSCName 
    reloadAndReinit [getMDBasename] false
    if {![isset binVelocities] && ![isset velocities]} {
        reinitvels [temperature]
    }
    return
}

# ::namdcph::initializeAlch
#
#   Initialize the alchemical settings needed for nonequilibrium alchemical
# switching of protonation states. This activates the appropriate nonbonded
# kernels and stores SimParameters settings that will be needed during 
# switches.
#
proc ::namdcph::initializeAlch {} {
    if {[isset alch]} {
        cphAbort "Constant pH is currently incompatible with alchemical"\
                 "transformations. Remove all alch settings."
    }
    cphPrint "Setting up nonequilibrium alchemical switching."
    # Here we require alchemical force computations and complete decoupling of
    # _all_ interactions. Unstaged linear coupling without softcore shifting is
    # currently assumed/enforced, but this might be ok to change.
    #
    alch on
    alchType TI
    alchDecouple off
    alchElecLambdaStart 0.0
    alchVdwLambdaEnd 1.0
    alchVdwShiftCoeff 0.0
    alchBondLambdaEnd 1.0
    alchBondDecouple on
    alchLambda 0.0
    alchLambda2 1.0
    alchOutFreq 0 ;# Suppress output - this would just be a nightmare.
    # Bonds between corresponding alchemical atoms are built dynamically when
    # a switch is started. In principle, this does not conflict with any other
    # extraBonds settings. However, indices in the extraBondsFile are likely
    # corrupted and would lead to wildly unexpected behavior - disable this for
    # now.
    #
    # TODO: Read and modify these files on they fly?
    if {[isset extraBonds]} {
        cphAbort "Constant pH is currently incompatible with extraBonds."
    }
    extraBonds on
    extraBondsFile [swFilename extrabonds]
    clearExtraBonds
    run 0
    # NB: Because alchLambda is incremented _before_ dynamics steps, if 
    # alchLambdaFreq were nonzero when "run 0" is called we would get a 
    # confusing energy log with a fractional alchLambda.
    #
    alchLambdaFreq 1
    alch off
    return
}

# ::namdcph::printSettingsSummary
#
# Print a summary of the constant pH settings.
#
proc ::namdcph::printSettingsSummary {} {
    variable ::namdcph::configFilename
    variable ::namdcph::restartFilename
    set StarBar "***************************************"
    cphPrint $StarBar
    cphPrint "CONSTANT pH MD ACTIVE"
    if {[string length $restartFilename] > 0} {
        cphPrint "RESTART FILENAME $restartFilename"
    }
    cphPrint "SYSTEM pH [getpH]"
    cphPrint "CONSTANT pH CONFIGURATION FILE $configFilename"
    cphPrint "NONEQUILIBRIUM SWITCH PARAMETERS:"
    cphPrint "ALCHEMICAL FORCE CONSTANT [getAlchFrcCons] kcal/mol-A^2"
    cphPrint "neMD/MC CRITERION TEMPERATURE [$::thermostatTempCmd]"
    cphPrint "TEMPORARY FILE INFORMATION:"
    cphPrint "cpH TOPOLOGY FILE BASENAME [getMDBasename]"
    cphPrint "neMD/MC TRAJECTORY BASENAME [getSWBasename]"
    cphPrint $StarBar
    return
}

# ::namdcph::printSystemSummary
#
# Print a summary of the titratable system. 
#
proc ::namdcph::printSystemSummary {} {
    set StarBar "***************************************"
    cphPrint $StarBar
    cphPrint "TITRATABLE RESIDUE DEFINITIONS:"
    cphPrint "[join [cphSystem get resdefs] " "]"
    cphPrint "TITRATABLE SYSTEM SUMMARY:"
    cphPrint "[cphSystem get numresidues] RESIDUE(S)"
    cphPrint "[llength [cphSystem get occupancy]] PROTONATION SITE(S)"
    cphPrint $StarBar
    cphPrint [format "%-19s : %5s : %s" segid:resid:resname state pKai]
    foreach resLabel [cphSystem get reslabel] state [cphSystem get state]\
            pKaiList [cphSystem get pKai] {
        cphPrint [format "%-19s : % 5s : %-s" $resLabel $state $pKaiList]
    }
    cphPrint $StarBar
    return
}

# ::namdcph::printTitratorSummary
#
# Print a summary of the MC moves.
#
proc ::namdcph::printTitratorSummary {} {
    set StarBar "*******************************************"
    cphPrint $StarBar
    cphPrint "CONSTANT pH neMD/MC MOVES:"
    cphPrint [format "%-19s : %8s %12s" "move label" numsteps weight]
    dict for {moveLabel data} [getMoveSet] {
        set numsteps [dict get $data numsteps]
        set weight [dict get $data weight]
        cphPrint [format "%-19s : % 8d % 12.2f" $moveLabel $numsteps $weight]
    }
    cphPrint $StarBar
    return
}

# ::namdcph::printProposalSummary
proc ::namdcph::printProposalSummary {segresidList} {
    set retList [list]
    foreach segresid $segresidList {
        set reslabel [cphSystem get reslabel $segresid]
        set state [cphSystem get state $segresid]
        set trialState [cphSystem get trialState $segresid]
        lappend retList "$reslabel:$state:$trialState"
#        cphPrint "Switch $reslabel $state --> $trialState"
    }
    return $retList
}

# ::namdcph::printTitratorReport
#
# Print a report of the titration MC statistics.
#
proc ::namdcph::printTitratorReport {} {
    set StarBar "*******************************************"
    cphPrint $StarBar
    cphPrint "CONSTANT pH MD STATISTICS:"
    cphPrint [format "%-19s : %-8s %-12s" "move label" attempts "accept. rate"]
    dict for {moveLabel data} [getMoveSet] {
        set attempts [dict get $data attempts]
        set successes [dict get $data successes]
        if {$successes && $attempts} {
            set rate [expr {100.*$successes/$attempts}]
        } else {
            set rate 0.0
        }
        cphPrint [format "%-19s : %8d %12.2f" $moveLabel $attempts $rate] 
    }
    cphPrint $StarBar
    return
}

# =============================================================================
# Getter Routines
#
# These are largely unnecessary, but cut down on "variable" declarations.
# =============================================================================
proc ::namdcph::getpH {} {
    return $::namdcph::SystempH
}

proc ::namdcph::getAlchFrcCons {} {
    return $::namdcph::AlchFrcCons
}

proc ::namdcph::getMDBasename {} {
    return $::namdcph::MDBasename
}

proc ::namdcph::getSWBasename {} {
    return $::namdcph::SWBasename
}
