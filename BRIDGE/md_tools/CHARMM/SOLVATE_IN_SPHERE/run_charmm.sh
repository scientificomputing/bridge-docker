#!/bin/sh

# _________ read inputs from the galaxy wrapper __________

crdin=$1 
psfin=$2 
radius=$3
s_origin=$4      
r_origin=$5
extrain=$6
top=$7
par=$8
crdout=$9
psfout=${10}
xplorpsfout=${11}
pdbout=${12}
wat=${13} 
output=${14}
tooldir=${15}

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

# _________create a stream for inputs_________  

streamfile=./variables.str
echo "*stream file for inputs" > $streamfile
echo "*                   " >> $streamfile
echo "                    " >> $streamfile
echo "set crdin \"$crdin\"" >> $streamfile
echo "set psfin \"$psfin\"" >> $streamfile
echo "set ffields \"$ffields\"" >> $streamfile
echo "set top \"$top\"" >> $streamfile
echo "set par \"$par\"" >> $streamfile
echo "set crdout \"$crdout\"" >> $streamfile
echo "set psfout \"$psfout\"" >> $streamfile
echo "set xplorpsfout \"$xplorpsfout\"" >> $streamfile
echo "set pdbout \"$pdbout\"" >> $streamfile
echo "                    " >> $streamfile
echo "return"   >> $streamfile

# _________copy or link all other inputs_________

cp ${tooldir}/toppar.str .

$charmm  << HEREDOC > $output
* FILENAME: solvate_in_sphere.inp
* PURPOSE:  solvate the system in a water (TIP3) water sphere 
* AUTHOR:   Tharindu Senapathi
* 
DIMENS CHSIZE 3000000 MAXRES 3000000   

BOMLEV -2        

!_____stream in variables_____    

stream variables.str

stream ./toppar.str


if $extrain .eq. yes  then

open read card unit 10 name @top
read  rtf card unit 10 append

open read card unit 20 name @par
read para flex card unit 20 append

endif


! Read in Topology and  Parameter files
set RADIUS $radius
set overlap 2.2

!SET LIMITS FOR FINDING CORRECT SPHERE  0,0,0 is the catalytic residue My protein is elongated along x-axis
set xstart -3 !We want to include more of the top of the protein
set ystart -3
set zstart -3  
set xlim 3 
set ylim 3 
set zlim 3
set sol all 
! default
if @?ADDWATERS ne 1 then
 set ADDWATERS 2 
endif


! origin
set orig segid $s_origin .and. resid $r_origin

echo @orig

open read card unit 14 name @psfin
read psf unit 14 card
close unit 14

open read card unit 13 name @crdin
read coord card unit 13
close unit 13

hbuild 

define origin sele @{orig} end

!____________ Move sphere around _________

coor orient sele all end !Rotate so that longest axis is on the x-axis
coor statistics sele origin end
coordinate translate xdirection -?XAVE ydirection -?YAVE zdirection -?ZAVE 

! Rename all waters to SPH0
rename segid SPH0 sele resname TIP3 end

define original sele (point 0.0 0.0 0.0 cut @RADIUS) end
echo ?nsel

!the next part searches for the best position of a potential sphere
!it first moves the sphere in the x direction to find the best position
!it then moves the sphere in the y direction to see if it can improve on the previous position
!then it moves the sphere in the z direction to see if it can improve on the previous position
 
! move sphere around
set prev 0
set x @xstart 
set maxx 0.0
label loop1
       define current sele ( (point @x 0.0 0.0 cut @RADIUS ) .and. .not. resname TIP3 ) end !select everything in current sphere
       set select ?nsel
       if @select gt @prev then !when the selection is greater than the previous one
       set prev @select
       set maxx @x
       else
       increment x by 1
       endif
       if @x le @xlim goto loop1
set y @ystart
set maxy 0.0
label loop2       
       define current sele ( (point @maxx @y 0.0 cut @RADIUS ) .and. .not. resname TIP3 ) end !select everything in current sphere
       set select ?nsel
       if @select gt @prev then !when the selection is greater than the previous one
       set prev @select
       set maxy @y
       else
       increment y by 1
       endif
if @y le @ylim goto loop2
set z @zstart
set maxz 0.0
label loop3
       define current sele ( (point @maxx @maxy @z cut @RADIUS ) .and. .not. resname TIP3 ) end !select everything in current sphere
       set select ?nsel
       if @select gt @prev then !when the selection is greater than the previous one
       set prev @select
       set maxz @z
       else
       increment z by 1
       endif
if @z le @zlim goto loop3
echo ?maxx ?maxy ?maxz
echo @prev
!Get the statistics 
coordinate statistics select origin end

!Center reaction center 

coor translate xdir -@maxx ydir -@maxy zdir -@maxz sele  all end ! Move the system 
set res ?nres


! delete the ions
delete atom sele segid SOD .or. segid CLA end

label minimize 

! define your solute
dele atom sele .byres. ( (.not. point 0.0 0.0 0.0 cut @RADIUS) .and. resname TIP3) end

MMFP  !This potential is similar to the sbound potential
GEO  sphere quartic -
     force 0.1 droff @RADIUS p1 2.25 select type OH2 end
END

cons fix sele .not. resname TIP3 end

delete atom sele segid SHP0 .or. resid TIP3 end

mini sd nstep 100
mini conj nstep 100

cons fix sele none end

open write unit 15 card name @crdout
write unit 15 coor card
close unit 15

open write unit 16 form name @psfout
write unit 16 psf card
close unit 16

open write unit 17 form name @pdbout
write coor pdb unit 17 
close unit 17

open write unit 18 form name @xplorpsfout
write psf xplor unit 18 card
close unit 18


stop

HEREDOC

cp  "${tooldir}/../../force_fields/sbmd_water_potentials/wat${radius}.pot" ${wat}
