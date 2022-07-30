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

set sphere_all [Sphere 900 $mid_plane_radius]
set sphere $sphere_all
puts "total no. of points [llength $sphere]"

set no_of_pts [llength $sphere]

proc check_cutoff {l cutoff} {
	set smaller_than_cutoff {}
	foreach element $l {
		if {$element < $cutoff} {
			lappend smaller_than_cutoff $element
		}
	}
	set len [llength $smaller_than_cutoff]
	return $len
} 

set pts_for_S {}
set i 0
set cutoff_S 100.0
set count_s 0
while {$i < 90} {
	set ind [expr int(rand()*$no_of_pts)]
	set pta [lindex $sphere $ind]
	set distance {}
	incr count_s
	if {[llength $pts_for_S] == 0} {
		lappend pts_for_S $pta
		set sphere [lreplace $sphere $ind $ind]
		set no_of_pts [llength $sphere]
		incr i
	} else {
		foreach pt $pts_for_S {
			set dist [vecdist $pt $pta]
			lappend distance $dist
		}
		if {[check_cutoff $distance $cutoff_S] == 0} {
			lappend pts_for_S $pta
	                set sphere [lreplace $sphere $ind $ind]
        	        set no_of_pts [llength $sphere]
 	               	incr i
		}
	}		
}



set pts_for_E {}
set j 0
set cutoff_E 80.0
set count_e 0
while {$j < 90} {
	set ind [expr int(rand()*$no_of_pts)]
	set pta [lindex $sphere $ind]
	set distance {}
	incr count_e
	if {[llength $pts_for_E] == 0} {
		lappend pts_for_E $pta
		set sphere [lreplace $sphere $ind $ind]
		set no_of_pts [llength $sphere]
		incr j
	} else {
		foreach pt $pts_for_E {
			set dist [vecdist $pt $pta]
			lappend distance $dist
		}
		if {[check_cutoff $distance $cutoff_E] == 0} {
			lappend pts_for_E $pta
	                set sphere [lreplace $sphere $ind $ind]
        	        set no_of_pts [llength $sphere]
 	               	incr j
		}
	}		
}


puts "count $count_s"
puts "no. of points for S [llength $pts_for_S]"

puts "count $count_e"
puts "no. of points for E [llength $pts_for_E]"

set pts_for_M $sphere
puts "no. of points for M [llength $pts_for_M]"


################################################################################################

set scent [vecscale [eval "vecadd $sphere_all"] [expr 1.0/[llength $sphere_all]]]
       
set count 0

proc prot_moving {prot_unit mid end pt scent} {
	set mid_cent [measure center $mid]
        set end_cent [measure center $end]
        set v1a [vecnorm [vecsub $end_cent $mid_cent]]
        set v1b [vecnorm [vecsub $pt $scent]]
        set ang1 [expr acos(max(-1.0,min(1.0,[vecdot $v1a $v1b])))]
        set cross1 "{[veccross $v1a $v1b]}"
	$prot_unit move [eval "trans center {$mid_cent} offset {$pt} axis $cross1 $ang1 rad"] 
}

for {set i 1} {$i <= 90 } {incr i} {
	set trimer [atomselect top "segname SA$i SB$i SC$i CA$i"]
        set end1 [atomselect top "segname SA$i SB$i SC$i and resid 1 to 305 and backbone"]
        set mid1 [atomselect top "segname SA$i SB$i SC$i and resid 1211 to 1236 and backbone"]
	set j [expr $i-1]
	set pt1 [lindex $pts_for_S $j]
        prot_moving $trimer $mid1 $end1 $pt1 $scent
	$trimer delete
	$mid1 delete
	$end1 delete
	incr count
	puts "$count done"	
}

puts "S protein done"

for {set i 1} {$i <= 90 } {incr i} {
	set pentamer [atomselect top "segname EA$i EB$i EC$i ED$i EE$i"]
        set end2 [atomselect top "segname EA$i EB$i EC$i ED$i EE$i and resid 35 to 40 and backbone"]
        set mid2 [atomselect top "segname EA$i EB$i EC$i ED$i EE$i and resid 18 to 38 and backbone"]
	set j [expr $i-1]
	set pt2 [lindex $pts_for_E $j]
        prot_moving $pentamer $mid2 $end2 $pt2 $scent
	$pentamer delete
	$mid2 delete
	$end2 delete
	incr count
	puts "$count done"	
}

puts "E protein done"

for {set i 1} {$i <= 720 } {incr i} {
	set dimer [atomselect top "segname L$i N$i"]
        set end3 [atomselect top "segname L$i N$i and resid 1 to 14 and backbone"]
        set mid3 [atomselect top "segname L$i N$i and ((resid 15 to 39) or (resid 75 to 99) or (resid 50 to 70)) and backbone"]
	set j [expr $i-1]
	set pt3 [lindex $pts_for_M $j]
        prot_moving $dimer $mid3 $end3 $pt3 $scent
	$dimer delete
	$mid3 delete
	$end3 delete
	incr count
	puts "$count done"	
}

puts "M protein done"

############### Saving the pdb and psf of liposome
[atomselect top all] writepdb S_E_M_dist_900d.pdb
[atomselect top all] writepsf S_E_M_dist_900d.psf


################### Finishing statement
puts "Finish At last"
##############################################################################################
