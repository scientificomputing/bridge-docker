#!/bin/sh

# _________ read inputs from the galaxy wrapper ____________

crdin=$1 
psfin=$2
buffer=$3
extrain=$4
top=$5
par=$6
setuppdb=$7
setupcrd=$8
setuppsf=$9
structureprm=${10}
output=${11}
tooldir=${12} 

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
echo "set buffer $buffer" >> $streamfile
echo "set ffields \"$ffields\"" >> $streamfile
echo "set top \"$top\"" >> $streamfile
echo "set par \"$par\"" >> $streamfile 
echo "set setuppdb \"$setuppdb\"" >> $streamfile 
echo "set setupcrd \"$setupcrd\"" >> $streamfile
echo "set setuppsf \"$setuppsf\"" >> $streamfile        
echo "set structureprm \"$structureprm\"" >> $streamfile
echo "                    " >> $streamfile
echo "return"   >> $streamfile

# _________copy or link all other inputs_________
# both  a link and a copy are dangerous , what if tools change. should be managed reference data --- FIXME . copy is also bad because duplicating. 
cp ${tooldir}/toppar.str . # both a link and a copy are dangerous , what if tools change. should be managed reference data
cp ${tooldir}/tip3.crd . 
cp ${tooldir}/crystal_image.str . 

$charmm  << HEREDOC > $output
* FILENAME: system_setup.inp
* PURPOSE:  setup the system for MD simulations.
* AUTHOR:   Adopted from CHARMM-GUI input scripts modified by Tharindu Senapathi
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


open unit 10 read form name @psfin
read psf unit 10 card

open unit 14 read card name @crdin
read coor unit 14 card

! determine the molecular extent
coor orient 
coor stat sele all end

calc Xinit = int ( ( ?Xmax - ?Xmin ) + 2 * @buffer ) + 1
calc Yinit = int ( ( ?Ymax - ?Ymin ) + 2 * @buffer ) + 1
calc Zinit = int ( ( ?Zmax - ?Zmin ) + 2 * @buffer ) + 1
calc Lbox = @Xinit
if Yinit .gt. @Lbox calc Lbox = @Yinit
if Zinit .gt. @Lbox calc Lbox = @Zinit

calc Xinit = @Lbox
calc Yinit = @Lbox
calc Zinit = @Lbox
calc Rinit = @Lbox / 2.0

calc xcen = 0.0
calc ycen = 0.0
calc zcen = 0.0

delete atom sele all end

! parameters for water box & image (crystal) 

calc BoxSizeX = @Xinit
calc BoxSizeY = @Yinit
calc BoxSizeZ = @Zinit 

set XTLtype = ORTHorhombic
if BoxSizeX .eq. @BoxSizeY if BoxSizeX .ne. @BoxSizeZ set XTLtype = TETRagonal
if BoxSizeX .eq. @BoxSizeY if BoxSizeX .eq. @BoxSizeZ set XTLtype = CUBIc
set A = @BoxSizeX 
set B = @BoxSizeY
set C = @BoxSizeZ 
set Alpha = 90.0
set Beta  = 90.0
set Gamma = 90.0

! a pre-equilibrated water cubic box with L=18.8560
set L 18.8560

! number of boxes along XYZ-directions
calc Xnum = int(@BoxSizeX/@L) + 1
calc Ynum = int(@BoxSizeY/@L) + 1
calc Znum = int(@BoxSizeZ/@L) + 1

! base unit of water box
read sequence TIP3 216
generate W000 setup noangle nodihedral

open read unit 10 card name tip3.crd
read coor unit 10 card 
close unit 10

coor stat sele type OH2 end
calc Lhalf = @L / 2.0
coor trans xdir 1.0 dist @Lhalf
coor trans ydir 1.0 dist @Lhalf
coor trans zdir 1.0 dist @Lhalf
coor stat sele type OH2 end

! planar water box unit (XY)

set J2  1
label DO_2
    set J1  1
    label DO_1

    calc wsegid = ( @J2 - 1 ) * @Xnum + @J1

        read sequence TIP3 216
        generate W@wsegid setup noangle nodihedral

    coor duplicate select segid W000 end select segid W@wsegid end

    calc X = @L * ( @J1 - 1 )  ! mult X by @J1
    calc Y = @L * ( @J2 - 1 )  ! mult Y by @J2

    coor trans xdir @X ydir @Y select segid W@wsegid end

    incr J1 by 1
    if J1 le @Xnum goto DO_1
incr J2 by 1
if J2 le @Ynum goto DO_2

delete atom sele .byres. ( ( type OH2 ) .and. -
                           ( prop X .gt. @BoxSizeX .or. -
                             prop Y .gt. @BoxSizeY ) ) end

define WAT sele .byres. type OH2 end
if ?nsel .eq. 0 stop ! ABNORMAL TERMINATION: Too small box size 

delete atom sele segid W000 end

open write unit 10 card name water_tmp.crd
write coor unit 10 card

delete atom sele all end

! generate water box by stacking planar water boxes along Z

set J3  1
label DO_3

    open read card unit 10 name water_tmp.crd
    read sequence coor card unit 10
    generate Wz@J3 setup warn noangle nodihedral 

    open read unit 10 card name water_tmp.crd
    read coor unit 10 card append

    calc Z = @L * ( @J3 - 1 )  ! mult Y by @J3
    coor trans zdir @Z select segid Wz@J3 end

incr J3 by 1
if J3 le @Znum goto DO_3

delete atom sele .byres. ( ( type OH2 ) .and. -
                           ( prop Z .gt. @BoxSizeZ ) ) end
                           
coor stat sele type OH2 end
coor orient norotation 
coor stat sele type OH2 end

!Shaping the bo

COOR CONVERT ALIGNED SYMMETRIC @A @B @C @alpha @beta @gamma
coor copy comp

CRYSTAL DEFINE @XTLtype @A @B @C @alpha @beta @gamma
CRYSTAL BUILD NOPER 0 CUTOFF 2.0

!Image centering by residue
IMAGE BYRESID XCEN @xcen YCEN @ycen ZCEN @zcen

nbond ctonnb 2.0 ctofnb 3.0 cutnb 3.0 cutim 3.0 wmin 0.001
CRYSTAL FREE

coor diff comp
delete atom sele .byres. ( ( prop Xcomp .ne. 0 ) .or. -
                           ( prop Ycomp .ne. 0 ) .or. -
                           ( prop Zcomp .ne. 0 ) ) end

COOR CONVERT SYMMETRIC ALIGNED @A @B @C @alpha @beta @gamma


coor stat sele type OH2 end
set nwater ?nsel

open write unit 10 card name waterbox.crd
write coor unit 10 card
* Equilibrated water 
*

open write card unit 92 name waterbox.pdb
write coor pdb  unit 92 
* Equilibrated water 
*

open write unit 90 card name waterbox.str
write title unit 90
* read sequence TIP3 @nwater
* generate WAT setup noangle nodihedral
*

open  write unit 90 card name @structureprm
write title unit 90
* set BoxType  = rect
* set XTLtype  = @XTLtype
* set A = @A
* set B = @B
* set C = @C
* set Alpha = @Alpha
* set Beta  = @Beta
* set Gamma = @Gamma
* set xcen = @xcen
* set ycen = @ycen
* set zcen = @zcen
*

delete atom sele all end

! _____ Step 2: Add ions to neutralize the system _____

! Read PSF and Coordinates
open read unit 10 card name @psfin
read psf  unit 10 card

open read unit 10 card name @crdin
read coor unit 10 card

! protein volume, calculation with a grid spacing of 0.5
coor orient 
coor stat

calc dcel = 0.5
calc xdim = int ( ( ?xmax - ?xmin + 5.0 ) / @dcel ) + 1
calc ydim = int ( ( ?ymax - ?ymin + 5.0 ) / @dcel ) + 1
calc zdim = int ( ( ?zmax - ?zmin + 5.0 ) / @dcel ) + 1
calc space = @xdim * @ydim * @zdim

scalar wmain = radius
scalar wmain add 1.4    ! use solvent accessible surface for molecular volume
scalar 1 = wmain 
scalar 2 set 6.0
coor volume hole space @space sele .not. resname TIP3 end

set molvol = ?volume

! system volume from parameters for water box
stream @structureprm

calc sysvol = @A * @B * @C


!
! conc   : concentration (M)
! volumn : ion accessible volume (Ang**3)
! npos   : number of positive ions
! nneg   : number of negative ions
!

calc conc   = 0.15 
calc volumn = @sysvol - @molvol
set posval   = 1
set negval   = 1

calc ncharge = int( ?cgtot )

calc npos = 0
calc nneg = 0
if ncharge .lt. 0 calc npos = int( abs( ?cgtot ) / @posval )
if ncharge .gt. 0 calc nneg = int( abs( ?cgtot ) / @negval )

calc npos = @npos + int ( @conc * 6.021 * 0.0001 * @volumn ) * @negval
calc nneg = @nneg + int ( @conc * 6.021 * 0.0001 * @volumn ) * @posval

calc diff   = int ( ?cgtot + @npos*@posval - @nneg*@negval )
label neutral
   if diff .gt. 0 calc nneg = @nneg + 1
   if diff .lt. 0 calc npos = @npos + 1
   calc diff   = int ( ?cgtot + @npos*@posval - @nneg*@negval )
if diff .ne. 0 goto neutral

if npos .lt. 0 stop ! something wrong
if nneg .lt. 0 stop ! something wrong

! Randomly place the ions

!Generate SOD
if npos .gt. 0 then
   read sequence SOD @npos
   generate SOD warn
endif

!Generate CLA
if nneg .gt. 0 then
   read sequence CLA @nneg
   generate CLA warn
endif

calc nion = @npos + @nneg
if nion .eq. 0 goto continue

cons fix sele .not. ( segid SOD .or. segid CLA ) end

! initial positions for all ions

calc xpos = @A / 2.0
calc ypos = @B / 2.0
calc zpos = @C / 2.0

coor set xdir @xpos ydir @ypos zdir @zpos select segid SOD .or. segid CLA end

! Initial placement of ions


calc i   = 1
label doinit

    if i .le. @nion then
       calc j = @i
    endif
    if i .gt. @nion then
       calc j = int( ?random * @nion ) + 1
    endif
    
    set ion = SOD
    if j .gt. @npos then
       set ion = CLA
       calc j  = @j - @npos
    endif
    
    define target select segid @ion .and. resid @j end
    coor stat sele target end
    calc xsave  = ?xave
    calc ysave  = ?yave
    calc zsave  = ?zave
    
    calc xpos = @A * ( ?random - 0.5 )
    calc ypos = @B * ( ?random - 0.5 )
    calc zpos = @C * ( ?random - 0.5 )
    
    ! check if the ions are too close to solute 
    coor set xdir @xpos  ydir @ypos  zdir @zpos select target end
    update
    
    coor dist cut 4.5 sele target end -
                      sele ( .not. target ) .and. .not. hydrogen end
    
    if ?npair .gt. 0 then
       coor set xdir @xsave  ydir @ysave  zdir @zsave select target end
       goto doinit
    endif

increase i by 1
if i .le. @nion goto doinit

!Images & Energy Setup

cons fix sele none end

COOR CONVERT ALIGNED SYMMETRIC @A @B @C @alpha @beta @gamma

open read unit 10 card name crystal_image.str
CRYSTAL DEFINE @XTLtype @A @B @C @alpha @beta @gamma
CRYSTAL READ UNIT 10 CARD

!Image centering by residue
IMAGE BYATOM XCEN @xcen YCEN @ycen ZCEN @zcen sele resname SOD .or. resname CLA end

!Note eps=80, truncation at 20 Angstrom
nbonds atom switch vatom vswitch -
       ctonnb 20.0 ctofnb 20.0 cutnb 21.0 cutim 21.0 -
       inbfrq -1 imgfrq -1 wmin 1.0 cdie eps 80.0

cons fix sele .not. ( segid SOD .or. segi CLA ) end

!
! Monte Carlo (MC) simulations of ions
! 

!Random number generation
RAND UNIF ISEED 1521208506

calc nmc = 2000

calc i   = 1
label domc 

    if i .le. @nion then
       calc j = @i
    endif       
    if i .gt. @nion then
       calc j = int( ?random * @nion ) + 1
    endif

    set ion = SOD
    if j .gt. @npos then
       set ion = CLA
       calc j  = @j - @npos
    endif

    define target select segid @ion .and. resid @j end 
    coor stat sele target end
    calc xsave  = ?xave
    calc ysave  = ?yave
    calc zsave  = ?zave

    calc xpos = @A * ( ?random - 0.5 )
    calc ypos = @B * ( ?random - 0.5 )
    calc zpos = @C * ( ?random - 0.5 )

    ! check if the ions are too close to solute
    coor set xdir @xpos  ydir @ypos  zdir @zpos select target end
    coor dist cut 4.5 sele target end -
                      sele ( .not. target ) .and. .not. hydrogen end
    if ?npair .gt. 0 then
       coor set xdir @xsave  ydir @ysave  zdir @zsave select target end
       goto domc
    endif

    ! before the move
    coor set xdir @xsave  ydir @ysave  zdir @zsave select target end

    interaction sele target end sele .not. target end
    set pener = 99999.0
    if ?ener .lt. @pener calc pener = ?ener

    ! after the move
    coor set xdir @xpos  ydir @ypos  zdir @zpos select target end
    update

    interaction sele target end sele .not. target end
    calc dener = ?ener - @pener
    calc boltz = exp ( - @dener / 0.59 )

    if xsave .eq. @a       goto next
    if dener .lt. 0.0      goto next
    if ?random .lt. @boltz goto next

    coor set xdir @xsave  ydir @ysave  zdir @zsave select target end
    update

    label next

incr i by 1
if i le @nmc goto domc

ENERGY
COOR CONVERT SYMMETRIC ALIGNED @A @B @C @alpha @beta @gamma
CRYSTAL FREE

open write unit 10 card name ions.pdb
write coor unit 10 pdb 
close unit 10

delete atom sele .not. ( segid SOD .or. segid CLA ) end

open write unit 10 card name ions.crd
write coor unit 10 card
close unit 10

label continue

open write unit 10 card name ions.prm
write title unit 10
* set npos = @npos ! Number of positive ions
* set nneg = @nneg ! Number of negative ions
* set posid = SOD
* set negid = CLA
*

open write unit 10 card name ions.str
write title unit 10
* read sequence SOD @npos !Generate SOD
* generate SOD warn
* read sequence CLA @nneg    !Generate CLA
* generate CLA warn
* 

delete atom sele all end

! _____ Step 3: Combine systems ______

! Read PSF and Coordinates
open unit 10 read form name @psfin
read psf unit 10 card

open unit 14 read card name @crdin
read coor unit 14 card

!Reorient Solute (should be here)
coor orient 
coor stat sele all end

!Read Water in
stream waterbox.str

open read card unit 30 name waterbox.crd
read coor card unit 30 append
close unit 30

! Add Ion

stream ions.prm

if npos .gt. 0 then
   read sequence @posid @npos
   generate @posid warn
endif

if nneg .gt. 0 then
   read sequence @negid @nneg
   generate @negid warn
endif

calc nion = @npos + @nneg

if nion .gt. 0 then
   open read unit 10 card name ions.crd
   read coor unit 10 card append
endif

!
! Remove water molecules close to or overlapped with 
! biomolecule, crystal water, and generated ions
!

define TOTO sele .not. hydrogen .and. .not. segid WAT end
coor stat sele TOTO end
coor stat sele type OH2 .and. segid WAT end

delete atom sele .byres. ( ( type OH2 .and. segid WAT ) .and. ( TOTO .around. 2.8 ) ) end

join WAT renumber
coor stat sele type OH2 .and. segid WAT end
set nwater ?nsel

mini abnr nstep 100

open unit 13 write card name @setupcrd
write coor card unit 13
* system setup
*

open write unit 14 card name @setuppsf
write psf  unit 14 card
* system setup
*

open write card unit 15 name @setuppdb
write coor pdb  unit 15 
* system setup
*

stop

stop

HEREDOC
