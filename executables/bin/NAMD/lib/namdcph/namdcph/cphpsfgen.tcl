# cphpsfgen.tcl
#
#   This file provides extensions to psfgen for handling alchemical PSFs, 
# specifically those related to protonation state transformations as are
# encountered in constant pH applications.
#
package require Tcl 8.5


variable BOLTZMANN 0.001987191 ;# in kcal/mol-K

if {[catch numPes]} {
    # numPes is a NAMD only command - this will catch for regular Tcl.
    proc print {args} {
        puts [join $args " "]
    }
    proc abort {args} {
        print $args
        exit
    }
}

if {![catch [package require psfgen 1.6.5]]} {
    abort "cphpsfgen requires psfgen 1.6.5 (shipped with NAMD 2.12) or later"
}
source [file join [file dirname [info script]] "numtcl.tcl"]


proc psfgenPrint {args} {
    print "cphpsfgen) [join $args " "]"
}

# =============================================================================
# Core Routines using psfgen
# =============================================================================
# alchPatch
#
#   Apply an alchemical patch. This assumes that the patch template has an even
# number of atoms divided into two groups (l0atoms and l1atoms) - the first is 
# fully coupled at alchLambda = 0 and the second is fully coupled at 
# alchLambda = 1. It is further assumed that there is a one-to-one ordered 
# mapping of the two groups. After the patch is applied, coordinates are copied
# from l0atoms to l1atoms. If desired, Gaussian noise is added under the 
# assumption of an isotropic harmonic potential between the atoms. Random 
# velocities for the new atoms are allocated based on the given temperaure and 
# their masses.
#
# NB: Unlike the standard patch command, this is currently limited to 
# application of a single patch.
#
# Arguments:
# ----------
# patch : string
#   Name of the patch
# segresid : string
#   Residue specification formatted as "<segid>:<resid>" - this is the same
#   syntax for the regular psfgen patch command.
# l0atomList : list
#   Names of atoms that are fully coupled at alchLambda = 0.
# l1atomList : list
#   Names of atoms that are fully coupled at alchLambda = 1.
# frcCons : float (optional, default 0.0)
#   Force constant to be used for adding Gaussian noise during copy placement 
#   of alchemical atoms. If a non-positive value is given, these will just be
#   copied without noise.
# temp : float (optional, default 0.0)
#   Temperature at which to allocate randomized coordinates and velocities. If
#   a non-positive value is given, these will just be copied.
# buildH : bool (optional, default false)
#   If true, do not copy hydrogens attached to carbon atoms; build them based 
#   on the IC table instead. This prevents catastrophic problems when enforcing
#   constraints on equivalent hydrogens.
#
# Returns:
# --------
# None
#
proc alchPatch {patch segresid l0atomList l1atomList 
                {frcCons 0.0} {temp 0.0} {buildH false}} {
    set kT [expr {$::BOLTZMANN*$temp}]
    set SigmaX [expr {($frcCons > 0) ? [expr {sqrt($kT/(2*$frcCons))}] : 0.0}]
    lassign [split $segresid ":"] segid resid
    set residue [segment residue $segid $resid]
    psfgenPrint "creating new alchemical side chain for $segresid:$residue"
    patch $patch $segresid
    # Copy/modify/reallocate atoms at lambda = 1 based on lambda = 0 atoms. 
    #
    set PrevElem ""
    set NSkipH 0
    foreach l0atom $l0atomList l1atom $l1atomList {
        set Elem [string index $l0atom 0]
        # Coordinates
        if {($buildH) && ($Elem == "H") 
            && (($PrevElem == "H") || ($PrevElem == "C") 
                || ($PrevElem == "N"))} {
            incr NSkipH
        } else {
            set x [segment coordinates $segid $resid $l0atom]
            AddNoiseToVector x 0.0 $SigmaX
            psfset coord $segid $resid $l1atom $x
        }
        set PrevElem $Elem
        # Velocities
        if {$temp > 0} {
            set Mass [segment mass $segid $resid $l0atom]
            set SigmaV [expr {($Mass > 0) ? [expr {sqrt($kT/$Mass)}] : 0.0}]
            set v {0.0 0.0 0.0}
            AddNoiseToVector v 0.0 $SigmaV
        } else {
            set v [segment velocities $segid $resid $l0atom]
        }
        psfset vel $segid $resid $l1atom $v
        # Alchemical Groups
        psfset beta $segid $resid $l0atom -1.0
        psfset beta $segid $resid $l1atom 1.0
    }
    if {($buildH) && ($NSkipH > 0)} {
        psfgenPrint "$NSkipH equivalent hydrogens not copied, building now..."
        guesscoord
    }
    return
}

# alchUnpatch
#
#   Remove one or more alchemical patches by merging the lambda = 1 group of
# atoms into the lambda = 0 group. The lambda = 0 group is _deleted_ and only 
# the _names_ of the second group are effected - coordinates and velocities are
# retained. Alchemical flags are reset to zero.
# 
# Arguments:
# ----------
# segresid : string (multiple arguments possible)
#   Residue specification formatted as "<segid>:<resid>" - this is the same
#   syntax for the regular psfgen patch command.
# l0atomList : list
#   Names of atoms that are fully coupled at alchLambda = 0.
# l1atomList : list
#   Names of atoms that are fully coupled at alchLambda = 1.
#
# Returns:
# --------
# None
#
proc alchUnpatch {segresid l0atomList l1atomList} {
    lassign [split $segresid ":"] segid resid
    set residue [segment residue $segid $resid]
    psfgenPrint "merging alchemical sidechains for $segresid:$residue"
    foreach l0atom $l0atomList l1atom $l1atomList { 
        psfset coord $segid $resid $l1atom\
                [segment coordinates $segid $resid $l0atom]
        delatom $segid $resid $l0atom
        psfset beta $segid $resid $l1atom 0.0
        psfset name $segid $resid $l1atom $l0atom
    }
    return
}

# psfgenRead
#
# Convenience routine for reading all possible psfgen input using a basename.
#
# NB: This also clears the current topology. 
# 
# Arguments:
# ----------
# basename : string
#   Basename of .psf, .pdb, .coor, and .vel files to read 
#
# Returns:
# --------
# None
#
proc psfgenRead {basename} {
    resetpsf
    readpsf "$basename.psf" pdb "$basename.pdb" namdbin "$basename.coor" \
            velnamdbin "$basename.vel"
    return
}

# psfgenWrite
#
# Convenience routine for writing all possible psfgen output using a basename.
# 
# Arguments:
# ----------
# basename : string
#   Basename of .psf, .pdb, .coor, and .vel files to write 
#
# Returns:
# --------
# None
#
proc psfgenWrite {basename {oldXSC ""}} {
    writepsf nopatches "$basename.psf"
    writepdb "$basename.pdb"
    writenamdbin "$basename.coor" velnamdbin "$basename.vel"
    if {[string length $oldXSC] > 0} {
        file copy -force $oldXSC "$basename.xsc"
    }
    return
}

#==============================================================================
# Math Routines
#==============================================================================
# AddNoiseToVector
#
# Add Gaussian noise to a vector with mean mu and width sigma.
#
# Arguments:
# ----------
# vector : list
#   The vector to be modified
# mu : float
#   mean of the noise
# sigma : float
#   standard deviation of the noise
#
# Returns:
# --------
# None
#
proc AddNoiseToVector {vector mu sigma} {
    upvar $vector v
    if {$sigma <= 0.0} {
        return
    }
    for {set k 0} {$k < [llength $v]} {incr k} {
        set w [lindex $v $k]
        set dw [expr {normal(0.0, $sigma)}]
        lset v $k [expr {$w + $dw}]
    }
    return
}

