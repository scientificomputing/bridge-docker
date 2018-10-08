# namdcph.tcl
#
# Main script for running constant pH in NAMD. 
#
#   This is essentially a dirty hack that tricks the user into renaming parts 
# of the constant pH namespace (defined in namdcph.tcl) into the global 
# namespace created by the NAMD interpreter. The advantage is that we get all
# of the tidiness and security of the private namespace without the user having
# even to be aware that it exists. All of the procedures declared in this file 
# are essentially interface functions that modify variables in another 
# namespace.
#
source [file join [file dirname [info script]] "namdcph/namdcph.core.tcl"]

namespace eval ::namdcphwrapper {
    namespace import ::namdcph::*
    namespace export cphRun runcph pH cphConfigFile cphOutFile\
            cphSetResidueState cphSetResiduepKai\
            cphRestartFile cphRestartFreq\
            cphForceConstant cphMDBasename cphSwitchBasename\
            cphMaxProposalAttempts cphNumMinSteps cphProposalWeight\
            cphNumstepsPerSwitch\
            testResidue cphExcludeResidue

    # Old command name for legacy reasons.
    proc runcph {args} {
        cphRun {*}$args
    }
}

namespace import ::namdcphwrapper::*
