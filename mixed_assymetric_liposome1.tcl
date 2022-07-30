
###########################################################################################
################# Initial input data

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

puts "$inner_num_vec"
puts "$outer_num_vec"

