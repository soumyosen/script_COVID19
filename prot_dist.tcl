############################################################################################
set inner_radius 367.0
set outer_radius 400.0
set PI 3.1415926535897931
################### Homogeneous distribution of points on the surface of a sphere

set mid_plane_radius [expr ($inner_radius + $outer_radius)/2.0]

proc Sphere {num_pts radius} {
         global PI
         set dlong [expr $PI*(3-sqrt(5.0))]
         set dz [expr 2.0/$num_pts]
         set long 0.0
         set z [expr 1-$dz/2]
         set sphere {}
         for {set k 1} {$k<=$num_pts} {incr k} {
                set r [expr sqrt(1-$z*$z)]
                set sp_x [expr $radius*cos($long)*$r]
                set sp_y [expr $radius*sin($long)*$r]
                set sp_z [expr $radius*$z]
                lappend sphere "$sp_x $sp_y $sp_z"
                set z [expr $z-$dz]
                set long [expr $long+$dlong]
         }
         return $sphere
}

set sphere_all [Sphere 600 $mid_plane_radius]
set sphere $sphere_all
puts "total no. of points [llength $sphere]"

set no_of_pts [llength $sphere]
set pts_for_S {}
for {set k 1} {$k <= 60} {incr k} {
	set ind [expr int(rand()*$no_of_pts)]
	set pta [lindex $sphere $ind]
	lappend pts_for_S $pta
	set sphere [lreplace $sphere $ind $ind]
	set no_of_pts [llength $sphere]
	#puts "$no_of_pts"
}

puts "no. of points for S [llength $pts_for_S]"


set pts_for_E {}
for {set k1 1} {$k1 <= 60} {incr k1} {
	set ind1 [expr int(rand()*$no_of_pts)]
	set ptb [lindex $sphere $ind1]
	lappend pts_for_E $ptb
	set sphere [lreplace $sphere $ind1 $ind1]
	set no_of_pts [llength $sphere]
}

puts "no. of points for E [llength $pts_for_E]"

set pts_for_M $sphere

puts "no. of points for M [llength $pts_for_M]"

################################################################################################

set scent [vecscale [eval "vecadd $sphere_all"] [expr 1.0/[llength $sphere_all]]]
       
set count 0
 
for {set i 1} {$i <= 60 } {incr i} {
	set trimer [atomselect top "segname SA$i SB$i SC$i CA$i"]
        set mid1 [atomselect top "segname SA$i SB$i SC$i and resid 1237 to 1273 and backbone"]
        set end1 [atomselect top "segname SA$i SB$i SC$i and resid 1 to 305 and backbone"]
        set memb_mid_pt1 [atomselect top "segname SA$i SB$i SC$i and resid 1211 to 1236 and backbone"]
	set mid1_cent [measure center $mid1]
	set end1_cent [measure center $end1]
	set memb_mid_pt1_cent [measure center $memb_mid_pt1]
	set j [expr $i-1]
	set pt1 [lindex $pts_for_S $j]
	set v1a [vecnorm [vecsub $end1_cent $mid1_cent]]
        set v1b [vecnorm [vecsub $pt1 $scent]]
        set ang1 [expr acos(max(-1.0,min(1.0,[vecdot $v1a $v1b])))]
        set cross1 "{[veccross $v1a $v1b]}"
        $trimer move [eval "trans center {$mid1_cent} axis $cross1 $ang1 rad"]
        $trimer move [eval "trans center {$memb_mid_pt1_cent} offset {$pt1}"]
        $trimer delete
	$mid1 delete
	$end1 delete
	$memb_mid_pt1 delete
	incr count
	puts "$count done"	
}

puts "S protein done"

for {set i 1} {$i <= 60 } {incr i} {
	set pentamer [atomselect top "segname EA$i EB$i EC$i ED$i EE$i"]
        set mid2 [atomselect top "segname EA$i EB$i EC$i ED$i EE$i and resid 1 to 17 and backbone"]
        set end2 [atomselect top "segname EA$i EB$i EC$i ED$i EE$i and resid 35 to 40 and backbone"]
        set memb_mid_pt2 [atomselect top "segname EA$i EB$i EC$i ED$i EE$i and resid 18 to 38 and backbone"]
	set mid2_cent [measure center $mid2]
	set end2_cent [measure center $end2]
	set memb_mid_pt2_cent [measure center $memb_mid_pt2]
	set j [expr $i-1]
	set pt2 [lindex $pts_for_E $j]
	set v1a [vecnorm [vecsub $end2_cent $mid2_cent]]
        set v1b [vecnorm [vecsub $pt2 $scent]]
        set ang1 [expr acos(max(-1.0,min(1.0,[vecdot $v1a $v1b])))]
        set cross1 "{[veccross $v1a $v1b]}"
        $pentamer move [eval "trans center {$mid2_cent} axis $cross1 $ang1 rad"]
        $pentamer move [eval "trans center {$memb_mid_pt2_cent} offset {$pt2}"]
        $pentamer delete
	$mid2 delete
	$end2 delete
	$memb_mid_pt2 delete
	incr count
	puts "$count done"	
}

puts "E protein done"

for {set i 1} {$i <= 480 } {incr i} {
	set dimer [atomselect top "segname L$i N$i"]
        set mid3 [atomselect top "segname L$i N$i and resid 201 to 222 and backbone"]
        set end3 [atomselect top "segname L$i N$i and resid 1 to 14 and backbone"]
        set memb_mid_pt3 [atomselect top "segname L$i N$i and ((resid 15 to 39) or (resid 75 to 99) or (resid 50 to 70)) and backbone"]
	set mid3_cent [measure center $mid3]
	set end3_cent [measure center $end3]
	set memb_mid_pt3_cent [measure center $memb_mid_pt3]
	set j [expr $i-1]
	set pt3 [lindex $pts_for_M $j]
	set v1a [vecnorm [vecsub $end3_cent $mid3_cent]]
        set v1b [vecnorm [vecsub $pt3 $scent]]
        set ang1 [expr acos(max(-1.0,min(1.0,[vecdot $v1a $v1b])))]
        set cross1 "{[veccross $v1a $v1b]}"
        $dimer move [eval "trans center {$mid3_cent} axis $cross1 $ang1 rad"]
        $dimer move [eval "trans center {$memb_mid_pt3_cent} offset {$pt3}"]
        $dimer delete
	$mid3 delete
	$end3 delete
	$memb_mid_pt3 delete
	incr count
	puts "$count done"	
}

puts "M protein done"

############### Saving the pdb and psf of liposome
[atomselect top all] writepdb S_E_M_dist_600.pdb
[atomselect top all] writepsf S_E_M_dist_600.psf


################### Finishing statement
puts "Finish At last"
##############################################################################################
