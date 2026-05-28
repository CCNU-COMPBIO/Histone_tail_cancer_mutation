# ---- User settings ----
set systems {WT}
set runs {run1 run2 run3}
set outdir "./DNA_contacts"

set cutoff 4.0
set time_step_ns 0.04
set traj_step 5
set eff_dt [expr {$time_step_ns * $traj_step}]

if {![file exists $outdir]} { file mkdir $outdir }

# ---- Main loop ----
foreach name $systems {
    foreach run $runs {


        if {![file exists $gro_file] || ![file exists $xtc_file]} {
            puts "WARNING: $gro_file or $xtc_file not found. Skipping."
            continue
        }

        puts "Processing system: $name, run: $run ..."
        mol delete all
        mol load gro $gro_file
        mol addfile $xtc_file step $traj_step waitfor all

        if {[file exists "rename_chainID_sym.tcl"]} {
            source rename_chainID_sym.tcl
        } else {
            puts "ERROR: rename_chainID_sym.tcl is missing."
            exit
        }

        set nframes [molinfo top get numframes]

        # Output files
        set f_all   [open "$outdir/dna_prot_contacts_${name}_${run}_all.dat" w]
        set f_h3    [open "$outdir/dna_prot_contacts_${name}_${run}_h3.dat" w]
        set f_h4    [open "$outdir/dna_prot_contacts_${name}_${run}_h4.dat" w]
        set f_h2a_n [open "$outdir/dna_prot_contacts_${name}_${run}_h2a_n.dat" w]
        set f_h2a_c [open "$outdir/dna_prot_contacts_${name}_${run}_h2a_c.dat" w]
        set f_h2b   [open "$outdir/dna_prot_contacts_${name}_${run}_h2b.dat" w]

        # Header
        foreach f [list $f_all $f_h3 $f_h4 $f_h2a_n $f_h2a_c $f_h2b] {
            puts -nonewline $f "Time"
            for {set r -73} {$r <= 72} {incr r} { puts -nonewline $f "\t$r" }
            puts -nonewline $f "\n"
        }

        for {set i 0} {$i < $nframes} {incr i} {
            set current_time [format "%.3f" [expr {$eff_dt * $i}]]
            if {[expr {$i % 100}] == 0} { puts "  Frame $i / $nframes" }

            foreach f [list $f_all $f_h3 $f_h4 $f_h2a_n $f_h2a_c $f_h2b] {
                puts -nonewline $f "$current_time"
            }

            # Histone selections (heavy atoms only)
            set sel_h3    [atomselect top "(segname CHA CHE) and (resid 1 to 36) and noh" frame $i]
            set sel_h4    [atomselect top "(segname CHB CHF) and (resid 1 to 20) and noh" frame $i]
            set sel_h2a_n [atomselect top "(segname CHC CHG) and (resid 1 to 13) and noh" frame $i]
            set sel_h2a_c [atomselect top "(segname CHC CHG) and (resid 119 to 129) and noh" frame $i]
            set sel_h2b   [atomselect top "(segname CHD CHH) and (resid 1 to 24) and noh" frame $i]

            # All-histone selection (matches the segments/ranges above)
            set sel_all [atomselect top "(((segname CHA CHE) and (resid 1 to 36)) or ((segname CHB CHF) and (resid 1 to 20)) or ((segname CHC CHG) and (resid 1 to 13 119 to 129)) or ((segname CHD CHH) and (resid 1 to 24))) and noh" frame $i]

            # DNA base-pair positions: r in [-73, 72], paired with -r in the opposite strand
            for {set r1 -73} {$r1 <= 72} {incr r1} {
                set r2 [expr {$r1 * (-1)}]
                set sel_dna [atomselect top "(((segname CHI) and (resid $r1)) or ((segname CHJ) and (resid $r2))) and noh" frame $i]

                set nc_all   [llength [lindex [measure contacts $cutoff $sel_all   $sel_dna] 1]]
                set nc_h3    [llength [lindex [measure contacts $cutoff $sel_h3    $sel_dna] 1]]
                set nc_h4    [llength [lindex [measure contacts $cutoff $sel_h4    $sel_dna] 1]]
                set nc_h2a_n [llength [lindex [measure contacts $cutoff $sel_h2a_n $sel_dna] 1]]
                set nc_h2a_c [llength [lindex [measure contacts $cutoff $sel_h2a_c $sel_dna] 1]]
                set nc_h2b   [llength [lindex [measure contacts $cutoff $sel_h2b   $sel_dna] 1]]

                puts -nonewline $f_all   "\t$nc_all"
                puts -nonewline $f_h3    "\t$nc_h3"
                puts -nonewline $f_h4    "\t$nc_h4"
                puts -nonewline $f_h2a_n "\t$nc_h2a_n"
                puts -nonewline $f_h2a_c "\t$nc_h2a_c"
                puts -nonewline $f_h2b   "\t$nc_h2b"

                $sel_dna delete
            }

            foreach f [list $f_all $f_h3 $f_h4 $f_h2a_n $f_h2a_c $f_h2b] {
                puts -nonewline $f "\n"
            }

            $sel_all delete
            $sel_h3 delete
            $sel_h4 delete
            $sel_h2a_n delete
            $sel_h2a_c delete
            $sel_h2b delete
        }

        foreach f [list $f_all $f_h3 $f_h4 $f_h2a_n $f_h2a_c $f_h2b] { close $f }
        puts "Done: $name, $run"
    }
}

puts "Finished."
exit
