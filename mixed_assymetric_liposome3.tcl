
set lipid {CHL1 POPC POPE POPS POPI POPA PSM POPI13 POPI24 POPI33}
set lipid_file {CHOL POPC POPE POPS POPI POPA PSM PIP1 PIP2 PIP3}
set outerlipid_ratio {34 39 6 0 0 0 21 0 0 0}
set innerlipid_ratio {29 17 26 11 3 2 9 1 1 1}
set lipid_SA {40.0 68.3 58.8 60.4 67.4 60.1 55.4 67.4 67.4 67.4}
set inner_radius 365
set outer_radius 400

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

set inner_num_vec [lipidNum $innerlipid_ratio $lipid_SA $inner_radius]
set outer_num_vec [lipidNum $outerlipid_ratio $lipid_SA $outer_radius]
puts "$lipid"
puts "$inner_num_vec"
puts "$outer_num_vec"

##############################################################################################
############# Make a list of residue numbers and shuffle the list randomly such that 
############# different lipids can be distributed on the spherical surface randomly

set resnum [llength [lsort -unique [[atomselect top all] get residue]]]
set resname_list [lsort -unique [[atomselect top all] get resname]]
puts "resnum $resnum"
puts "resname list $resname_list"

proc list_shuffle {l} {
	set num [llength $l]
	for {set i 0} {$i < $num} {incr i} {
		set j [expr {int(rand()*$num)}]
		set temp [lindex $l $j]
		set l [lreplace $l $j $j [lindex $l $i]]
		set l [lreplace $l $i $i $temp]
	}
	return $l
}

set residue_name {}
set residue_count {}
set residue_num_list {}
set inner_count {}
set inner_num_list {}

foreach residue_type $resname_list {
	set num  [lsort -unique [[atomselect top "resname $residue_type"] get residue]]
	set count [llength $num]
	lappend residue_name $residue_type
	lappend residue_count $count
	for {set k 0} {$k < $count} {incr k} {
		lappend residue_num_list [lindex $num $k]
	}
	set ind [lsearch $lipid $residue_type]
	set inner_layer_lipid_num [lindex $inner_num_vec $ind] 
	lappend inner_count $inner_layer_lipid_num
	for {set s 0} {$s < $inner_layer_lipid_num} {incr s} {
		lappend inner_num_list [lindex $num $s]
	} 
}

puts "residue name $residue_name"
puts "residue count $residue_count"
puts "residue inner count $inner_count"

set outer_num_list {}
foreach num_any $residue_num_list {
	if { [ lsearch $inner_num_list $num_any ] < 0 } {lappend outer_num_list $num_any}
}

puts "inner count [llength $inner_num_list]"
puts "outer count [llength $outer_num_list]"

set fp1 [open "trial_inner_num.dat" w]
set fp2 [open "trial_inner_num_shuffled.dat" w]
set fp3 [open "trial_outer_num.dat" w]
set fp4 [open "trial_outer_num_shuffled.dat" w]

set inner_num_list_shuffled [list_shuffle $inner_num_list]
set outer_num_list_shuffled [list_shuffle $outer_num_list]
#puts "$inner_num_list_shuffled"

puts $fp1 "$inner_num_list"
puts $fp2 "$inner_num_list_shuffled"
puts $fp3 "$outer_num_list"
puts $fp4 "$outer_num_list_shuffled"
