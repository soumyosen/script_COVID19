
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

set inner_num_vec [lipidNum $lipid_ratio $lipid_SA $inner_radius]
set outer_num_vec [lipidNum $lipid_ratio $lipid_SA $outer_radius]

puts "$inner_num_vec"
puts "$outer_num_vec"

############################################################################################
############### Provide name of the lipid and number of lipids in inner and outer sphere
############### At the end give the combined pdb and psf of all the lipids with 
############### whatever numbers are calculated by previous function

proc genpsfpdb {a in out} {
        mol new ../$a.psf
        mol addfile ../$a.pdb
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

foreach lpf $lipid_file in_num $inner_num_vec out_num $outer_num_vec {
        genpsfpdb $lpf $in_num $out_num
}

writepsf combined.psf
writepdb combined.pdb
exec rm -f a*
