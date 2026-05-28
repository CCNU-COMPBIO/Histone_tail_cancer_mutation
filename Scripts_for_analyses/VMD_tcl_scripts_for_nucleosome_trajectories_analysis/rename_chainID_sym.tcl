proc renumber { sel start } {
	if { [$sel num] == 0 } {
		puts "Error in renumber: empty selection!"
		return
	}
	set oresid [ $sel get resid ]
	set delta [ expr $start - [ lindex $oresid 0] ]
	set nresid { }
	foreach r $oresid {
		lappend nresid [ expr $r + $delta ]
	}
	$sel set resid $nresid
}


# 选择并修改链A的残基
set sel [atomselect top "resid 1 to 146 "]
$sel set chain I
$sel set segname CHI
renumber $sel -73
$sel delete

# 选择并修改链B的残基
set sel [atomselect top "resid 147 to 292 "]
$sel set chain J
$sel set segname CHJ
renumber $sel -72
$sel delete

# 选择并修改链C的残基
set sel [atomselect top "resid 293 to 427 "]
$sel set chain A
$sel set segname CHA
renumber $sel 1
$sel delete

# 选择并修改链D的残基
set sel [atomselect top "resid 428 to 562 "]
$sel set chain E
$sel set segname CHE
renumber $sel 1
$sel delete

# 选择并修改链E的残基
set sel [atomselect top "resid 563 to 664 "]
$sel set chain B
$sel set segname CHB
renumber $sel 1
$sel delete

# 选择并修改链F的残基
set sel [atomselect top "resid 665 to 766 "]
$sel set chain F
$sel set segname CHF
renumber $sel 1
$sel delete

# 选择并修改链G的残基
set sel [atomselect top "resid 767 to 895 "]
$sel set chain C
$sel set segname CHC
renumber $sel 1
$sel delete

# 选择并修改链H的残基
set sel [atomselect top "resid 896 to 1024 "]
$sel set chain G
$sel set segname CHG
renumber $sel 1
$sel delete

# 选择并修改链I的残基
set sel [atomselect top "resid 1025 to 1149 "]
$sel set chain D
$sel set segname CHD
renumber $sel 1
$sel delete

# 选择并修改链J的残基
set sel [atomselect top "resid 1150 to 1274 "]
$sel set chain H
$sel set segname CHH
renumber $sel 1
$sel delete

