# ---- User settings ----
set systems {WT}
set runs {run1 run2 run3}
set cutoff 4.0

set time_step_ns 0.04
set traj_step 5
set eff_dt [expr {$time_step_ns * $traj_step}]

set outdir "./tail_contacts"
if {![file exists $outdir]} { file mkdir $outdir }

# ---- Main loop ----
foreach name $systems {
    foreach run $runs {
        puts "Processing system: $name, run: $run ..."
        mol delete all

        if {![file exists $gro_file] || ![file exists $xtc_file]} {
            puts "WARNING: $gro_file or $xtc_file not found. Skipping."
            continue
        }

        mol load gro $gro_file
        mol addfile $xtc_file step $traj_step waitfor all

        if {[file exists "rename_chainID_sym.tcl"]} {
            source rename_chainID_sym.tcl
        } else {
            puts "ERROR: rename_chainID_sym.tcl is missing."
            exit
        }

        set nframes [molinfo top get numframes]

        # ---- Output helper ----
        proc open_and_header {filepath resid_start resid_end resid_extra_start resid_extra_end} {
            set out [open $filepath w]
            puts -nonewline $out "Time"
            for {set r $resid_start} {$r <= $resid_end} {incr r} { puts -nonewline $out "\t$r" }
            if {$resid_extra_start != ""} {
                for {set r $resid_extra_start} {$r <= $resid_extra_end} {incr r} { puts -nonewline $out "\t$r" }
            }
            puts -nonewline $out "\n"
            return $out
        }

        # ---- Open output files ----
        # H3 (1-36)
        set f11 [open_and_header "$outdir/tail_contacts_all_${name}_${run}_h3_A.dat" 1 36 "" ""]
        set f12 [open_and_header "$outdir/tail_contacts_all_${name}_${run}_h3_E.dat" 1 36 "" ""]
        set f13 [open_and_header "$outdir/tail_contacts_backbone_${name}_${run}_h3_A.dat" 1 36 "" ""]
        set f14 [open_and_header "$outdir/tail_contacts_backbone_${name}_${run}_h3_E.dat" 1 36 "" ""]
        set f15 [open_and_header "$outdir/tail_contacts_base_${name}_${run}_h3_A.dat" 1 36 "" ""]
        set f16 [open_and_header "$outdir/tail_contacts_base_${name}_${run}_h3_E.dat" 1 36 "" ""]

        # H4 (1-20)
        set f21 [open_and_header "$outdir/tail_contacts_all_${name}_${run}_h4_B.dat" 1 20 "" ""]
        set f22 [open_and_header "$outdir/tail_contacts_all_${name}_${run}_h4_F.dat" 1 20 "" ""]
        set f23 [open_and_header "$outdir/tail_contacts_backbone_${name}_${run}_h4_B.dat" 1 20 "" ""]
        set f24 [open_and_header "$outdir/tail_contacts_backbone_${name}_${run}_h4_F.dat" 1 20 "" ""]
        set f25 [open_and_header "$outdir/tail_contacts_base_${name}_${run}_h4_B.dat" 1 20 "" ""]
        set f26 [open_and_header "$outdir/tail_contacts_base_${name}_${run}_h4_F.dat" 1 20 "" ""]

        # H2A (1-13, 119-129)
        set f31 [open_and_header "$outdir/tail_contacts_all_${name}_${run}_h2a_C.dat" 1 13 119 129]
        set f32 [open_and_header "$outdir/tail_contacts_all_${name}_${run}_h2a_G.dat" 1 13 119 129]
        set f33 [open_and_header "$outdir/tail_contacts_backbone_${name}_${run}_h2a_C.dat" 1 13 119 129]
        set f34 [open_and_header "$outdir/tail_contacts_backbone_${name}_${run}_h2a_G.dat" 1 13 119 129]
        set f35 [open_and_header "$outdir/tail_contacts_base_${name}_${run}_h2a_C.dat" 1 13 119 129]
        set f36 [open_and_header "$outdir/tail_contacts_base_${name}_${run}_h2a_G.dat" 1 13 119 129]

        # H2B (1-24)
        set f41 [open_and_header "$outdir/tail_contacts_all_${name}_${run}_h2b_D.dat" 1 24 "" ""]
        set f42 [open_and_header "$outdir/tail_contacts_all_${name}_${run}_h2b_H.dat" 1 24 "" ""]
        set f43 [open_and_header "$outdir/tail_contacts_backbone_${name}_${run}_h2b_D.dat" 1 24 "" ""]
        set f44 [open_and_header "$outdir/tail_contacts_backbone_${name}_${run}_h2b_H.dat" 1 24 "" ""]
        set f45 [open_and_header "$outdir/tail_contacts_base_${name}_${run}_h2b_D.dat" 1 24 "" ""]
        set f46 [open_and_header "$outdir/tail_contacts_base_${name}_${run}_h2b_H.dat" 1 24 "" ""]

        set all_files [list $f11 $f12 $f13 $f14 $f15 $f16 \
                            $f21 $f22 $f23 $f24 $f25 $f26 \
                            $f31 $f32 $f33 $f34 $f35 $f36 \
                            $f41 $f42 $f43 $f44 $f45 $f46]

        # ---- Frame loop ----
        for {set i 0} {$i < $nframes} {incr i} {
            set current_time [format "%.3f" [expr {$eff_dt * $i}]]
            if {[expr {$i % 100}] == 0} {
                puts "  Frame $i / $nframes (t = $current_time ns)"
            }

            foreach f $all_files { puts -nonewline $f "$current_time" }

            # DNA selections (heavy atoms only)
            set sel_dna_all  [atomselect top "nucleic and noh" frame $i]
            set sel_dna_bb   [atomselect top "nucleic and backbone and noh" frame $i]
            set sel_dna_base [atomselect top "nucleic and not backbone and noh" frame $i]

            # H3 (CHA / CHE), resid 1-36
            for {set r 1} {$r <= 36} {incr r} {
                set sel_A [atomselect top "(segname CHA and resid $r) and noh" frame $i]
                set sel_E [atomselect top "(segname CHE and resid $r) and noh" frame $i]

                puts -nonewline $f11 "\t[llength [lindex [measure contacts $cutoff $sel_dna_all  $sel_A] 1]]"
                puts -nonewline $f13 "\t[llength [lindex [measure contacts $cutoff $sel_dna_bb   $sel_A] 1]]"
                puts -nonewline $f15 "\t[llength [lindex [measure contacts $cutoff $sel_dna_base $sel_A] 1]]"

                puts -nonewline $f12 "\t[llength [lindex [measure contacts $cutoff $sel_dna_all  $sel_E] 1]]"
                puts -nonewline $f14 "\t[llength [lindex [measure contacts $cutoff $sel_dna_bb   $sel_E] 1]]"
                puts -nonewline $f16 "\t[llength [lindex [measure contacts $cutoff $sel_dna_base $sel_E] 1]]"

                $sel_A delete
                $sel_E delete
            }

            # H4 (CHB / CHF), resid 1-20
            for {set r 1} {$r <= 20} {incr r} {
                set sel_B [atomselect top "(segname CHB and resid $r) and noh" frame $i]
                set sel_F [atomselect top "(segname CHF and resid $r) and noh" frame $i]

                puts -nonewline $f21 "\t[llength [lindex [measure contacts $cutoff $sel_dna_all  $sel_B] 1]]"
                puts -nonewline $f23 "\t[llength [lindex [measure contacts $cutoff $sel_dna_bb   $sel_B] 1]]"
                puts -nonewline $f25 "\t[llength [lindex [measure contacts $cutoff $sel_dna_base $sel_B] 1]]"

                puts -nonewline $f22 "\t[llength [lindex [measure contacts $cutoff $sel_dna_all  $sel_F] 1]]"
                puts -nonewline $f24 "\t[llength [lindex [measure contacts $cutoff $sel_dna_bb   $sel_F] 1]]"
                puts -nonewline $f26 "\t[llength [lindex [measure contacts $cutoff $sel_dna_base $sel_F] 1]]"

                $sel_B delete
                $sel_F delete
            }

            # H2A (CHC / CHG), resid 1-13 and 119-129
            foreach r {1 2 3 4 5 6 7 8 9 10 11 12 13 119 120 121 122 123 124 125 126 127 128 129} {
                set sel_C [atomselect top "(segname CHC and resid $r) and noh" frame $i]
                set sel_G [atomselect top "(segname CHG and resid $r) and noh" frame $i]

                puts -nonewline $f31 "\t[llength [lindex [measure contacts $cutoff $sel_dna_all  $sel_C] 1]]"
                puts -nonewline $f33 "\t[llength [lindex [measure contacts $cutoff $sel_dna_bb   $sel_C] 1]]"
                puts -nonewline $f35 "\t[llength [lindex [measure contacts $cutoff $sel_dna_base $sel_C] 1]]"

                puts -nonewline $f32 "\t[llength [lindex [measure contacts $cutoff $sel_dna_all  $sel_G] 1]]"
                puts -nonewline $f34 "\t[llength [lindex [measure contacts $cutoff $sel_dna_bb   $sel_G] 1]]"
                puts -nonewline $f36 "\t[llength [lindex [measure contacts $cutoff $sel_dna_base $sel_G] 1]]"

                $sel_C delete
                $sel_G delete
            }

            # H2B (CHD / CHH), resid 1-24
            for {set r 1} {$r <= 24} {incr r} {
                set sel_D [atomselect top "(segname CHD and resid $r) and noh" frame $i]
                set sel_H [atomselect top "(segname CHH and resid $r) and noh" frame $i]

                puts -nonewline $f41 "\t[llength [lindex [measure contacts $cutoff $sel_dna_all  $sel_D] 1]]"
                puts -nonewline $f43 "\t[llength [lindex [measure contacts $cutoff $sel_dna_bb   $sel_D] 1]]"
                puts -nonewline $f45 "\t[llength [lindex [measure contacts $cutoff $sel_dna_base $sel_D] 1]]"

                puts -nonewline $f42 "\t[llength [lindex [measure contacts $cutoff $sel_dna_all  $sel_H] 1]]"
                puts -nonewline $f44 "\t[llength [lindex [measure contacts $cutoff $sel_dna_bb   $sel_H] 1]]"
                puts -nonewline $f46 "\t[llength [lindex [measure contacts $cutoff $sel_dna_base $sel_H] 1]]"

                $sel_D delete
                $sel_H delete
            }

            foreach f $all_files { puts -nonewline $f "\n" }

            $sel_dna_all delete
            $sel_dna_bb delete
            $sel_dna_base delete
        }

        foreach f $all_files { close $f }
        puts "Done: $name, $run"
    }
}

puts "Finished."
exit
