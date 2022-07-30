
###########################################################################################
################# Initial input data

set lipid {DOPC DOPE DOPS SSM CHL}
set lipid_ratio {108 9 36 27 20}
set lipid_SA {69.7 63.4 71.6 55.3 40.0}
set inner_radius 154
set outer_radius 190

#############################################################################################
################ Number of different lipid molecules on a specific spherical surface 
################ Output: n dimensional vector: where n is number of different lipid molecules

proc lipidNum {ratio SA radius} {
        set tot 0.0
        foreach r $ratio sa $SA {
                set rsa [expr $r*$sa]
                set tot [expr $tot+$rsa]
        }

        global PI
        set PI 3.1415926535897931
        set totarea [expr 4.0*$PI*pow($radius,2)]
        set unit [expr $totarea/$tot]

        set number {}
        foreach r $ratio {
                set each_num [expr int(floor([expr $r*$unit]))]
                lappend number $each_num
        }
        return $number
}

set inner_num_vec [lipidNum $lipid_ratio $lipid_SA $inner_radius]
set outer_num_vec [lipidNum $lipid_ratio $lipid_SA $outer_radius]

puts "$inner_num_vec"
puts "$outer_num_vec"

############################################################################################
############### Provide name of the lipid and number of lipids in inner and outer sphere
############### At the end give the combined pdb and psf of all the lipids with 
############### whatever numbers are calculated by previous function

proc genpsfpdb {a in out} {
        mol new $a.psf
        mol addfile $a.pdb
        set sel [atomselect top all]
        set tot_lipid [expr $in+$out]
        if {$tot_lipid >= 10000} {
                set mult [expr int([expr {ceil($tot_lipid/10000.0)}])]
                set t 0
                set alp [split "$a" {}]
                set seg [join [list [lindex $alp 2] [lindex $alp 3]] ""]
                for {set i 1} {$i <= $mult} {incr i} {
                        set begin [expr 10000*($i-1)]
                        set end [expr (10000*$i)-1]
                        for {set j $begin} {$j <= [expr min($end,[expr $tot_lipid-1])]} {incr j} {
                                $sel set segname $seg$i
                                set k [expr $j-$begin]
                                $sel set resid $k
                                set t [expr $t+1]
                                $sel writepdb a$a$t.pdb
                                $sel writepsf a$a$t.psf
                        }
                }
        } else {
                for {set i 1} {$i <= $tot_lipid} {incr i} {
                        $sel set segname $a
                        $sel set resid $i
                        $sel writepdb a$a$i.pdb
                        $sel writepsf a$a$i.psf
                }
        }
        $sel delete
        package require psfgen
        for {set i 1} {$i <= $tot_lipid} {incr i} {
                readpsf a$a$i.psf
                coordpdb a$a$i.pdb
        }

        return
}

foreach lp $lipid in_num $inner_num_vec out_num $outer_num_vec {
        genpsfpdb $lp $in_num $out_num
}

writepsf combined.psf
writepdb combined.pdb
exec rm -f a*
mol delete all
mol new combined.psf
mol addfile combined.pdb

################################################################################################
############### Calculation of total number of lipid molecules in the inner sphere and outer 
############### sphere

proc numpoints {ptvec} {
	set num 0
	foreach a $ptvec {
		set num [expr $num+$a]
	}
	return $num
}

set inner_num [numpoints $inner_num_vec]
set outer_num [numpoints $outer_num_vec]

puts "$inner_num"
puts "$outer_num"

###########################################################################################
############### Getting the coordinates on a spherical surface with a given radius with number 
############### of points which should be equal to the number of lipid molecules

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

set inner_sphere [Sphere $inner_num $inner_radius]
set outer_sphere [Sphere $outer_num $outer_radius]

##########################################################################################
############## Functions to move the lipid molecules to a specific position of 
############## inner sphere and outer sphere

proc innermoving {res pt cent} {
        set residue [atomselect top "residue $res"]
        set mid [measure center [atomselect top "residue $res and (type PL or type OHL)"]]
        set end [measure center [atomselect top "residue $res and type CTL3"]]
        set v1 [vecnorm [vecsub $end $mid]]
        set v2 [vecnorm [vecsub $pt $cent]]
        set ang [expr acos(max(-1.0,min(1.0,[vecdot $v1 $v2])))]
        set cross "{[veccross $v1 $v2]}"
        $residue move [eval "trans center {$mid} offset {$pt} axis $cross $ang rad"]
        $residue delete
        return
}

proc outermoving {res pt cent} {
        set residue [atomselect top "residue $res"]
        set mid [measure center [atomselect top "residue $res and (type PL or type OHL)"]]
        set end [measure center [atomselect top "residue $res and type CTL3"]]
        set v1 [vecnorm [vecsub $mid $end]]
        set v2 [vecnorm [vecsub $pt $cent]]
        set ang [expr acos(max(-1.0,min(1.0,[vecdot $v1 $v2])))]
        set cross "{[veccross $v1 $v2]}"
        $residue move [eval "trans center {$mid} offset {$pt} axis $cross $ang rad"]
        $residue delete
        return
}

##############################################################################################
############# Calculating the centers of inner and outer sphere

set inner_cent [vecscale [eval "vecadd $inner_sphere"] [expr 1.0/[llength $inner_sphere]]]
set outer_cent [vecscale [eval "vecadd $outer_sphere"] [expr 1.0/[llength $outer_sphere]]]

##############################################################################################
############# Make a list of residue numbers and shuffle the list randomly such that 
############# different lipids can be distributed on the spherical surface randomly

set resnum [llength [lsort -unique [[atomselect top all] get residue]]]

set numlist {}
for {set i 0} {$i < $resnum} {incr i} {
	lappend numlist $i
}

for {set i 0} {$i < $resnum} {incr i} {
	set j [expr {int(rand()*$resnum)}]
	set temp [lindex $numlist $j]
	set numlist [lreplace $numlist $j $j [lindex $numlist $i]]
	set numlist [lreplace $numlist $i $i $temp]
}

#################################################################################################
############ Move the lipid molecules

set count 1
for {set i 0} {$i < $resnum} {incr i} {
	if {$count <= $inner_num} {
		set pt [lindex $inner_sphere $i]
		set spc_res [lindex $numlist $i]
		innermoving $spc_res $pt $inner_cent
		set count [expr $count+1]
	} else {
		set pt [lindex $outer_sphere [expr $i-$inner_num]]
		set spc_res [lindex $numlist $i]
		outermoving $spc_res $pt $outer_cent
		set count [expr $count+1]
	}
}

#################################################################################################

[atomselect top all] writepdb mixed_liposome.pdb
[atomselect top all] writepsf mixed_liposome.psf

#################################################################################################

mol delete all
mol new mixed_liposome.psf
mol addfile mixed_liposome.pdb
[atomselect top "segname SSM and type OIL"] set type OHL
[atomselect top all] writepdb mixed_liposome.pdb
[atomselect top all] writepsf mixed_liposome.psf
puts "Finish At Last"



