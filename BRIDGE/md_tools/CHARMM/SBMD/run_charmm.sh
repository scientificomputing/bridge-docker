#!/bin/sh

# _________ read inputs from the galaxy wrapper ____________ 

crdin=$1
psfin=$2
wat=$3
radius=$4
buffer=$5
temp=$6
pro=$7
sub1=$8
sub2=$9
sub3=${10}
restart=${11}
rstin=${12}
extrain=${13}
top=${14}
par=${15} 
seed=${16}
dcd_freq=${17}
freq=$(($dcd_freq * 1000))
simulation_time=${18}
steps=$(($simulation_time * 1000))
region=${19}
pdbout=${20}
crdout=${21}
rstout=${22}
dcdout=${23}
output=${24}  
tooldir=${25}

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
echo "set wat \"$wat\"" >> $streamfile
echo "set radius \"$radius\"" >> $streamfile
echo "set buffer \"$buffer\"" >> $streamfile
echo "set ffields \"$ffields\"" >> $streamfile  
echo "set rstin \"$rstin\"" >> $streamfile
echo "set top \"$top\"" >> $streamfile  
echo "set par \"$par\"" >> $streamfile 
echo "set region \"$region\"" >> $streamfile 
echo "set pdbout \"$pdbout\"" >> $streamfile 
echo "set crdout \"$crdout\"" >> $streamfile 
echo "set rstout \"$rstout\"" >> $streamfile 
echo "set dcdout \"$dcdout\"" >> $streamfile
echo "                    " >> $streamfile
echo "return"   >> $streamfile 

# _________copy or link all other inputs_________  

cp ${tooldir}/toppar.str . 

$charmm  << HEREDOC > $output
* FILENAME: SBMD.inp
* PURPOSE:  Molecular dynamics using SBC
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
close unit 20

open read unit 19 form name @crdin
read coor card unit 19
close unit 19

define prot sele (segid $pro ) end 
define sub sele (segid $sub1 .or. segid $sub2 .or. segid $sub3 ) end
define water sele resname TIP3 end
define backbone sele prot .and. (type CA .or. type C .or. type N .or. type O .or. type OT*) end

! LANGEVIN REGIONS

calc buffer = int ( $radius - $buffer ) 

echo @buffer

set buf  @buffer        !BUFFER BOUNDARY  IF YOU CHANGE THIS YOU NEED TO REGENERATE THE SCREENING 
set lang $radius	!LANGEVIN BOUNDARY   IF YOU CHANGE THIS YOU NEED TO REGENERATE THE SCREENING 
set mid  @buffer        ! Extra: extend reaction region for protein but limits it if too far into buffer

! START OF SETTING UP BUFFER REGION FOR LANGEVIN DYNAMICS 

!Set columns in comparison set to 0.0

scalar xcomp set 0.0
scalar ycomp set 0.0
scalar zcomp set 0.0

! SET REACTION CENTER

! @buf and @mid is used only side chains  can stretch into buffer region, but not further than @mid

! DEFINITIONS OF REGIONS

define reaction sele ( .byres. ( point 0.0 0.0 0.0 cut @buf ) ) .and. .not.  ( backbone .and. .not. - 
point 0.0 0.0 0.0 cut @mid)  END  !reaction is done by residue

define reservoir sele (.not. point 0.0 0.0 0.0 cut @lang) .and. .not. water   end 

define buffer sele .not. (reaction .or. reservoir)  end
define OUT sele (reservoir .and. .bonded. buffer) end
define IR1 sele (buffer .and. .bonded. OUT)  end
define IR2 sele (buffer .and. .bonded. IR1)  end
define FIX13 select .bygroup. (IR1 .or. IR2)  end

scalar xcomp set 1.0 sele reaction end 
scalar xcomp store 1 !flag in x of comparison set
write coor pdb name regionfile.reaction.pdb sele recall 1 end

! SET BUFFER REGION 

scalar ycomp set 1.0 sele buffer .and. .not. (FIX13)  end 
scalar ycomp store 2


! PROTEIN AND LIGAND IN BUFFER
scalar zcomp set 1.0 sele recall 2 .and. .not. ( hydrogen .or.  water) end

scalar zcomp store 3

! Reservoir region

scalar wcomp set 1.0 sele .not. ( recall 1 .or. recall 2 ) end
scalar wcomp store 4

open unit 13 write form name @region
write coor card comp unit 13
* regions
* column 1: reaction region  by residue partioning 
* column 2: Buffer region atoms, any atoms (byres) 
* column 3: Protein Langenvin atoms (same as col. 2 but no H or tip3).
* column 4: Outer region atoms 
*

! ____________ step 2 : dynamics setup ______________

!CONTROL DYNAMICS

set temp   $temp   ! In K 

!Miscellaneous variables
set potential normal    
!set cutpot  30 !hopefully this is fine
!Variables for definitions


! FORCE CONSTANTS

set bfnonox 1.30        !STOCHASTIC BOUNDARY FORCE X CONSTANT(1) FOR NON-OXYGEN BACKBONE ATOMS
set bfox 1.22           !STOCHASTIC BOUNDARY FORCE X CONSTANT(3) FOR OXYGEN BACKBONE ATOMS
set bfside 0.73         !STOCHASTIC BOUNDARY FORCE X CONSTANT(2) FOR SIDE CHAIN ATOMS

! FRICTION CONSTANTS 

set fprot 200.0 ! protein
set fwat  62.0  ! friction constant water

! Region file name
set reg @region

cons fix sele none end

open unit 3 name @reg read form
read coor comp card unit 3
close unit 3

scalar xcomp store 1
scalar ycomp store 2
scalar zcomp store 3
scalar wcomp store 4

! 1:Reaction region containing everything in @buf 
! 2:buffer region
! 3:Every atom in buffer without water and hydrogens We do not couple langevin to hydrogens
! 4:Reservoir
 

! SET STOCHASTIC BOUNDARY FORCES 

scalar wmain set 0.0 sele all end 

! assign force for protein atoms
! _________________________________
!
! potential
!
!          f= 3 R T*< delta xi**2 >
!
! and < delta xi**2 > comes from
!
!                      8 pi**2
! < delta xi**2 > = ---------------
!                        3 B
!
! __________________________________
!

BOMLEV -1

scalar wmain set 1.0 sele recall 3 end  !VERY IMPORTANT

! ALL OTHER BACKBONE ATOM Force CONSTANTS
scalar wmain mult @bfnonox sele prot .AND. ( recall 3 .and. (type C -
.or. type CA .or. type N ) ) end 

! BACKBONE OXYGEN Force CONSTANT
scalar wmain mult @bfox sele prot .AND. ( recall 3 .and.  (type O .or. type OT*))  end

! SIDE CHAIN Force CONSTANT
scalar wmain mult @bfside sele prot .AND. ( recall 3 .and. .not. - 
( type C .or. type CA .or. type N .or. type O .or. type OT*) )  end !These factors are very general


! WE NEED TO MULTIPLY A FEW THINGS
! THE EQUATION SHOULD BE 4*(pi **2)kT/B=force
! Biophys J. 2009 April 22; 96(8): 3074–3081. 8*(pi**2)kt/B

scalar wmain mult 0.001987191 sele recall 3 end !boltzmann constant for a mole
scalar wmain mult @temp sele recall 3 end 
scalar wmain mult 3. sele all end 
scalar wmain store 5

! setting up of friction coefficients 

scalar xcomp set 0.0 sele all end
scalar ycomp set 0.0 sele all end
scalar zcomp set 0.0 sele all end
scalar wcomp set 0.0 sele all end

scalar xcomp set 1.0 sele recall 3 end  !Everything except water and hydrogens in buffer
scalar zcomp set 1.0 sele recall 3 .or. type OH2 end 

scalar ycomp recall 5 
scalar zcomp mult @fprot sele prot .AND. recall 3 end 

! WATER FRICTION CONSTANTS
scalar zcomp mult @fwat sele resname TIP3 .and. type OH2 end 
scalar wcomp recall 1 

scalar xcomp store 8 

! SCREEN FORCES AND FRICTION COEFFICENTS
! At the top this is calculated as in this formula
! The force constant is scaled by a screening function so that the atoms at the boundary feel a stronger restoring force than further from the boundary
! ____________________________________________

! Screening potential in buffer region
!
!         (r-rb)**2 * (3*re -rb - 2*r)
!  s(r) = ----------------------------
!                2 * (re-rb)**3
!  rb = 24 A
!  re = 30 A
!
!  r = rb  => s(r)=0.0
!  r = re  => s(r)=0.5
! _____________________________________________
!

! setup screening

set rbinc 0.5
set r @buf
set pr 0

label setscreening
calc s = ((@r-@buf)**2 * (3 * @lang - @buf -2 * @r)) / (2 * (@lang-@buf)**3)

scalar ycomp mult @s sele ( ( recall 8 .and. -
(point 0.0 0.0 0.0 cut @r ) .and. .not. -
( point 0.0 0.0 0.0 cut @pr ) .or. resname -
TIP3P ) )  end

scalar zcomp mult @s sele ( ( recall 8 .and. -
(point 0.0 0.0 0.0 cut @r ) .and. .not. -
( point 0.0 0.0 0.0 cut @pr ) .or. resname -
TIP3P ) )  end

calc pr = @r
increment r by @rbinc

if @r .LE. @lang GOTO setscreening


! PUT A GENERALIZED RESTRAINING FORCE ON THE PROTEIN AND LIGAND IN BUFFER 
 
cons harmonic force 1.5 exponent 2 sele recall 8 .and. backbone end 
 
! We leave the side chain to move 

cons fix sele recall 4 - 
.or. (.not. .byres. ( point 0.0 0.0 0.0 cut @lang ) ) .and. .not. resname TIP3 end !FIX EVERYTHING IN RESERVOIR

scalar ycomp store 2 !transfer values 2 is force constants 3 is frictions constants
scalar zcomp store 3

! SET FINAL FORCES AND FRICTION COEFFICIENTS 

scalar CONST recall 2  !THIS IS THE CONSTRAINT FORCES 
scalar FBETA recall 3

open read unit 4 formatted name @wat 
sbound read unit 4
close unit 4

sbound set xref 0.0 yref 0.0 zref 0.0 assign 1 -
  sele   (resname TIP3 .and. type OH2) end

! _________ Step 3 : MOLECULAR DYNAMICS START _________________

update  GROUP  SWITCH CDIE  VDW VSWI  GRAD QUAD -
CUTNB 15.0  CTOFNB 14.0 CTONNB 10.0  WMIN 1.5  EPS 1.0  

faster 1        
SHAKE BONH PARAM 

if $restart .eq. no then    

! small minimization 
mini abnr nstep 100 tolegr 0.0001

open unit 22 write unform name @dcdout 
open unit 20 write form name   @rstout
open write card unit 24 name   nvt.eng     
 
DYNA VERL strt nstep  $steps timest 0.001 rbuf @buf iseed  $seed-   
         firstt @temp finalt $temp tbath $temp twindh 5.0 twindl -5.0  tstruct $temp - !Temperature 
         iasvel 1 iasors 1 ISVFRQ 10 ieqfrq 10 ICHECW 1 -
         IHTFRQ 0 TEMINC 0 -
         ilbfrq 5 inbfrq -1 - !inbfrq found to be approx 10  we use ntrfrq 100 to remove center of mass motion
!File manager
         IUNWRI 20 IUNCRD 22 KUNIT 24 IUNVEL -1 -
         nprint $freq $freq 100 nsavc $freq
 
endif

if $restart .eq. yes then  

open unit 22 read form name @rstin

open unit 24 write unform name @dcdout
open unit 20 write form name @rstout
open unit 28 write  form name  nvt.eng

DYNA VERL restrt nstep  $steps timest 0.001 rbuf @buf iseed $seed -   
         firstt $temp finalt $temp tbath $temp twindh 5.0 twindl -5.0  tstruct $temp - !Temperature 
         nprint $freq iprfrq 100 nsavc $freq -
         iasvel 1 iasors 1 ISVFRQ 10 ieqfrq 10 ICHECW 1 -
         IHTFRQ 0 TEMINC 0 -
         ilbfrq 5 inbfrq -1 - !inbfrq found to be approx 10  we use ntrfrq 100 to remove center of mass motion
!File manager
         IUNREAD 22 IUNWRI 20 IUNCRD 24 KUNIT 28 IUNVEL -1 

endif

open write card unit 14 name @crdout 
write coor card unit 14
close unit 14

open write form unit 15 name @pdbout
write coor pdb unit 15
close unit 15 

STOP

EOF

HEREDOC
