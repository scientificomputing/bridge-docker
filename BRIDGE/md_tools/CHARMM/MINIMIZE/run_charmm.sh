#!/bin/sh

# _________ read inputs from the galaxy wrapper ____________ 

crdin=$1
psfin=$2
structureprm=$3
steps=$4
ref=$5
pro=$6
sub1=$7
sub2=$8
sub3=$9
extrain=${10}
top=${11}
par=${12}
pdbout=${13}
crdout=${14}
psfout=${15}
xplorpsfout=${16}
pmespec=${17}
bbrmsd=${18}
scrmsd=${19}
subsrmsd=${20}
output=${21}
tooldir=${22}

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
echo "set ffields \"$ffields\"" >> $streamfile
echo "set extrain \"$extrain\"" >> $streamfile
echo "set top \"$top\"" >> $streamfile
echo "set par \"$par\"" >> $streamfile 
echo "set pdbout \"$pdbout\"" >> $streamfile  
echo "set crdout \"$crdout\"" >> $streamfile
echo "set psfout \"$psfout\"" >> $streamfile
echo "set xplorpsfout \"$xplorpsfout\"" >> $streamfile            
echo "set pmespec \"$pmespec\"" >> $streamfile            
echo "set bbrmsd \"$bbrmsd\"" >> $streamfile            
echo "set scrmsd \"$scrmsd\"" >> $streamfile            
echo "set subsrmsd \"$subsrmsd\"" >> $streamfile                                              
echo "                    " >> $streamfile
echo "return"   >> $streamfile


# _________copy or link all other inputs_________

cp ${tooldir}/toppar.str . 
cp ${tooldir}/crystal_image.str . 
cp ${tooldir}/generatefft.py . 

$charmm  << HEREDOC > $output
* FILENAME: minimize.inp
* PURPOSE:  Setup periodic boundary conditions and energy minimization
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

DELETE ATOM SELE .BYRES. (SEGID WAT .AND. TYPE OH2 .AND. -
       ( ( .NOT. SEGID WAT .AND. .NOT. (HYDROGEN .OR. LONE ) ) -
       .AROUND. 2.5 ) ) END !delete waters closer than 2.5A from system

! waterbox parameters

stream @structureprm

COOR CONVERT ALIGNED SYMMETRIC @A @B @C @alpha @beta @gamma

open read unit 10 card name crystal_image.str
CRYSTAL DEFINE @XTLtype @A @A @A @alpha @beta @gamma
CRYSTAL READ UNIT 10 CARD

!Image centering by residue

IMAGE BYRESID XCEN @xcen YCEN @ycen ZCEN @zcen sele resname TIP3 end
IMAGE BYRESID XCEN @xcen YCEN @ycen ZCEN @zcen sele ( segid SOD .or. segid CLA ) end

! PME parameters

system "python generatefft.py ?XTLA ?XTLB ?XTLC > $pmespec"

stream @pmespec

!Setting all the cutoffs and Ewald for nonbonded calculations

nbonds atom vatom vfswitch bycb -
ctonnb 10.0 ctofnb 12.0 cutnb 16.0 cutim 16.0 -
inbfrq -1 imgfrq -1 wmin 1.0 cdie eps 1.0 -
ewald pmew fftx @fftx ffty @ffty fftz @fftz  kappa 0.32 spline order 6

energy

mini sd nstep $steps
mini abnr nstep $steps

open write unit 15 form name @pdbout
write coor pdb unit 15

open write unit 18 card name @crdout
write coor card unit 18

open write unit 16 form name @psfout
write psf unit 16 card

open write unit 17 form name @xplorpsfout
write psf xplor unit 17 card

if $ref .eq. yes  then 
define BB   sele segid $pro .and. (type C .or. type CA .or. type N .or. type O .or. type OT*) end
define SC   sele .not. BB .and. .not. hydrogen .and. ( segid $pro ) end
define SUBS  sele segid $sub1 .or. segid $sub2 .or. segid $sub3 end

scalar wmain set 0 sele all end
scalar wmain set 1 sele BB end
write coor pdb name @bbrmsd

scalar wmain set 0 sele all end
scalar wmain set 1 sele SC end
write coor pdb name @scrmsd

scalar wmain set 0 sele all end
scalar wmain set 1 sele SUBS end
write coor pdb name @subsrmsd
endif

stop

HEREDOC
