#!/bin/sh

# _________ read inputs from the galaxy wrapper ____________ 

crdin=$1
psfin=$2
structureprm=$3
pmeshspec=$4
temp=$5
pressure=$6
restart=$7
rstin=$8
restraints=$9
pro=${10}
sub1=${11}
sub2=${12}
sub3=${13}
bb=${14}
sc=${15}
subs=${16}
extrain=${17}
top=${18}
par=${19}
seed=${20}
dcd_freq=${21}
freq=$(($dcd_freq * 1000))
simulation_time=${22}
steps=$(($simulation_time * 1000))
dcdoutput=${23}
pdbout=${24}
crdout=${25}
rstout=${26}
dcdout=${27}
output=${28}
tooldir=${29}

charmm="${tooldir}/../../bin/CHARMM/charmm"

if [ -e $charmm ]
then
    echo INFO: Found a valid CHARMM binary
else
    echo ERROR: failed to find CHARMM binary
    echo INFO: Request a license at http://charmm.chemistry.harvard.edu/charmm_lite.php
    exit 1
fi

ln -s "${tooldir}/../../force_fields/toppar" ffields
ffields=./ffields


# _________create a stream file for inputs_________  

streamfile=./variables.str
echo "*stream file for inputs" > $streamfile
echo "*                   " >> $streamfile
echo "                    " >> $streamfile
echo "set crdin \"$crdin\"" >> $streamfile
echo "set psfin \"$psfin\"" >> $streamfile
echo "set structureprm \"$structureprm\"" >> $streamfile
echo "set pmeshspec \"$pmeshspec\"" >> $streamfile
echo "set ffields \"$ffields\"" >> $streamfile
echo "set rstin \"$rstin\"" >> $streamfile
echo "set extrain \"$extrain\"" >> $streamfile
echo "set top \"$top\"" >> $streamfile
echo "set par \"$par\"" >> $streamfile
echo "set pdbout \"$pdbout\"" >> $streamfile 
echo "set crdout \"$crdout\"" >> $streamfile 
echo "set rstout \"$rstout\"" >> $streamfile 
echo "set dcdout \"$dcdout\"" >> $streamfile
echo "                    " >> $streamfile
echo "return"   >> $streamfile
 
# _________copy or link all other inputs_________

cp ${tooldir}/toppar.str . 

$charmm  << HEREDOC > $output
* FILENAME: NPT_MD.inp
* PURPOSE:  Molecular dynamics in NPT ensemble
* AUTHOR:   Tharindu Senapathi
* 

DIMENS CHSIZE 3000000 MAXRES 3000000

BOMLEV -2
! _______ stream in variables_______

stream variables.str

stream ./toppar.str


if $extrain .eq. yes  then

open read card unit 10 name @top
read  rtf card unit 10 append

open read card unit 20 name @par
read para flex card unit 20 append

endif


open read unit 20 card name @psfin
read psf card unit 20

open read unit 19 form name @crdin
read coor card unit 19


! waterbox parameters

stream @structureprm


! crystal build for PBC, building all transformations with a minimum atom-atom contact distance of less than 20.00 Angstroms.
CRYSTAL DEFINE @XTLtype @A @A @A @alpha @beta @gamma
crystal build cutoff 20.0 noper 0
IMAGE BYRESID XCEN @xcen YCEN @ycen ZCEN @zcen sele resname TIP3 end
IMAGE BYRESID XCEN @xcen YCEN @ycen ZCEN @zcen sele ( segid SOD .or. segid CLA ) end

! PME parameters

stream @pmeshspec

!Setting all the cutoffs and Ewald for nonbonded calculations


shake bonh param fast

mini SD   nstep 20 nprint 5
mini ABNR nstep 20 nprint 5

! estimate Pmass from SYSmass (total system mass)
! [there could be problems with exreme values, such as  Pmass << SYSmass or Pmass >> SYSmass
scalar mass stat
calc Pmass = int ( ?stot  /  50.0 )


if $restraints .eq. yes  then

define BB   sele segid $pro .and. (type C .or. type CA .or. type N .or. type O .or. type OT*) end
define SC   sele .not. BB .and. .not. hydrogen .and. ( segid $pro ) end
define SUBS  sele segid $sub1 .or. segid $sub2 .or. segid $sub3 end

cons harm force $bb sele BB end
cons harm force $sc sele SC end
cons harm force $subs sele SUBS end

endif

set temp = $temp

open unit 20 write form name @rstout

if $dcdoutput .eq. TRUE 
open unit 22 write unform name @dcdout
set dcdunit = 22
endif  

if $dcdoutput .eq. FALSE 
set dcdunit = -1
endif

if $restart .eq. yes then
set runopt = rest
open unit 7 read form name  @rstin
set rstunit = 7
endif

if $restart .eq. no then
set runopt = strt
set rstunit = -1
endif


DYNA VERL @runopt iseed $seed -
TIMEstep 0.001 NSTEp $steps -
INBFrq 10 IMGFrq 10 IHTFREQ 0 IPRFrq 100 IXTFrq 1000 -
IEQFrq 100 NTRFrq 100 NSAVC $freq NPRINT $freq  ISVFRQ 100 -  !NTRFRQ 1000
cpt hoover -
pcons pref $pressure pmass @Pmass  pgamma 20 -
reft @TEMP firstt @TEMP finalt @TEMP tstruct @TEMP twindh 5.0 twindl -5.0 -
IASORS 1 IASVEL 1 -
elec atom ewald -
kappa 0.32 PMEwald spline order 6 fftx @fftx ffty @ffty fftz @fftz qcor 0.0 -
fswitch cdie - 
vdw vatom VFswitch BYGR - 
cutnb 14.0 ctonnb 10.0 ctofnb 12.0 eps 1.0 -
wmin 1.5 - 
NBSCale 0.8 - 
nbxmod 5 e14fac 1.0 - 
IUNREAD @rstunit IUNWRI 20 IUNCRD @dcdunit KUNIT -1 IUNVEL -1


open write unit 10 card name @pdbout
write coor unit 10 pdb

open write unit 10 card name @crdout
write coor unit 10 card
close unit 10

stop

HEREDOC

