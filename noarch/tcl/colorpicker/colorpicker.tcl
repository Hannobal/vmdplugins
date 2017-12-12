
package provide colorpicker 1.0

global env

namespace eval ::ColorPicker:: {
  variable w               ;# window handle
  variable colors {}
  variable paletterows 28
  variable squaresize 25
  variable selected 0
  variable rectwdth [expr $squaresize*3]
  variable recthght 30 ;#[expr 2*$squaresize]
  variable texthght 20
  variable selectorrows [expr int((1.0*$paletterows*$squaresize)/$rectwdth)]
  variable showlabels 0
  variable colorsort [list \
    blue blue2 blue3 iceblue cyan cyan2 cyan3 green green2 green3 lime \
    yellow yellow2 yellow3 orange orange2 orange3 red red2 red3 purple pink \
    magenta magenta2 mauve violet violet2 white gray silver black tan ochre \
  ]
}

proc ::ColorPicker::rgb2vmd {r g b} {
  set or [expr $r/255.0]
  set og [expr $g/255.0]
  set ob [expr $b/255.0]
  return [list $or $og $ob]
}

proc ::ColorPicker::vmd2hex {r g b} {
  set r [expr round(255.0*$r)]
  set g [expr round(255.0*$g)]
  set b [expr round(255.0*$b)]
  return [format #%02x%02x%02x $r $g $b]
}

proc ::ColorPicker::rgb2hex {r g b} {
  return [format #%02x%02x%02x $r $g $b]
}

proc ::ColorPicker::rgb2hsv {r g b} {
  set min [expr min($r,$g,$b)]
  set max [expr max($r,$g,$b)]
  set v $max
  if {$v==0} {
    return [list 0 0 0]
  }
  set s [expr 255*($max-$min)/$v]
  if {$s==0} {
    return [list 0 $s $v]
  }
  if {$max == $r} {
    return [list [expr 43*($g-$b)/($max-$min)] $s $v]
  } elseif {$max == $g} {
    return [list [expr 85+43*($b-$r)/($max-$min)] $s $v]
  } else {
    return [list [expr 171+43*($r-$g)/($max-$min)] $s $v]
  }
}

proc ::ColorPicker::readcolors {filename {append 0}} {
  variable colors
  set fp [open $filename r]
  if {$append==0} { set colors {} }
  while { [gets $fp line] >= 0 } {
      if {[regexp {^\s*\S+\s+\d+\s+\d+\s+\d+} $line]} {
      set lst [regexp -all -inline -- {\S+} $line]
      lassign $lst nm r g b
      set lst [concat $lst [::ColorPicker::rgb2hsv $r $g $b]]
      lappend colors $lst
    }
  }
  close $fp
  return
}

proc ::ColorPicker::colorpicker {} {
  global env
  variable w
  variable colors
  variable colorsort
  variable paletterows
  variable squaresize
  variable selected
  variable rectwdth
  variable recthght
  variable texthght
  variable selectorrows
  variable showlabels
  
  readcolors [file join $env(VMDDIR) plugins noarch tcl colorpicker colors]
  
  if { [winfo exists ".colorpicker"] } {
    wm deiconify $w
    return
  }
  
  set w [toplevel ".colorpicker"]
  wm title $w "Color Picker"
  wm resizable $w 0 0 
  
  # color selector
  set mx [llength $colorsort]
  canvas $w.imgselect -width [expr $rectwdth*$selectorrows] -height [expr ($texthght+$recthght)*(1+ceil($mx/$selectorrows))]
  set i 0
  set l 0
  set vis "hidden"
  if {$showlabels} {
    set vis "normal"
  }
  foreach col $colorsort {
    set idx [colorinfo index $col]
    lassign [colorinfo rgb $col] r g b
    set txtcol #000000
    if {[expr max($r,$g,$b)]<0.5} {set txtcol #ffffff}
    set l [expr $i/$selectorrows]
    set p [expr $i%$selectorrows]
    $w.imgselect create rectangle [expr $rectwdth*$p] [expr $l*($recthght+$texthght)] [expr $rectwdth*($p+1)] [expr ($l+1)*$recthght+$l*$texthght] -fill [::ColorPicker::vmd2hex $r $g $b] -width 1 -tag r$i
    $w.imgselect create text [expr int($rectwdth*(0.5+$p))] [expr $l*$texthght+$recthght*(1+$l)] -text "$idx $col" -font "Helvetica, 8" -tag l$i -anchor n
    $w.imgselect create text [expr int($rectwdth*(0.5+$p))] [expr $l*$texthght+$recthght*(0.5+$l)] -text "" -font "Helvetica, 8" -tag t$i -state $vis -fill $txtcol
    incr i
  }
  set l [expr $selected/$selectorrows]
  set p [expr $selected%$selectorrows]
  $w.imgselect create rectangle [expr $rectwdth*$p] [expr ($texthght+$recthght)*$l] [expr $rectwdth*($p+1)] [expr ($texthght+$recthght)*($l+1)] -width 4 -tags marker
  pack $w.imgselect -side top -fill x -anchor w
  bind $w.imgselect <Button-1> { ::ColorPicker::selectcolor %x %y }
  
  frame $w.frm
  button $w.frm.opn -text "Open palette" -command ::ColorPicker::changepalette
  grid $w.frm.opn -column 0 -row 0
  button $w.frm.add -text "Add palette" -command ::ColorPicker::addpalette
  grid $w.frm.add -column 1 -row 0
  button $w.frm.sort -text "Sort palette" -command ::ColorPicker::sortpalette
  grid $w.frm.sort -column 2 -row 0
  button $w.frm.updt -text "Update" -command ::ColorPicker::updateselector
  grid $w.frm.updt -column 3 -row 0
  button $w.frm.prnt -text "Save colors" -command ::ColorPicker::printcolorscript
  grid $w.frm.prnt -column 4 -row 0
  checkbutton $w.frm.dsplbl -indicatoron yes -text "Display color labels" -compound none -width 0 -border 0 -variable ::ColorPicker::showlabels -command ::ColorPicker::setlabels
  grid $w.frm.dsplbl -column 5 -row 0
  pack $w.frm -side top -fill x -anchor w
  
  # color palette
  set i 0
  set l 0
  frame $w.frm2
  canvas $w.frm2.imgpalette -width [expr $squaresize*$paletterows] -yscrollcommand "$w.frm2.scroll set"
  redrawpalette
  bind $w.frm2.imgpalette <Button-1> { ::ColorPicker::setcolor %x %y }
  bind $w.frm2.imgpalette <Button-4> { %W yview scroll -1 units }
  bind $w.frm2.imgpalette <Button-5> { %W yview scroll 1 units }
  grid $w.frm2.imgpalette -column 0 -row 0 -sticky nswe
  scrollbar $w.frm2.scroll -command "$w.frm2.imgpalette yview"
  grid $w.frm2.scroll -column 1 -row 0 -sticky ns
  pack $w.frm2 -side top -fill x -anchor w
}

proc ::ColorPicker::updateselector {} {
  variable w
  variable colorsort
  variable colors
  set i 0
  foreach col $colorsort {
    set idx [colorinfo index $col]
    lassign [colorinfo rgb $col] r g b
    set new [::ColorPicker::vmd2hex $r $g $b]
    set old [$w.imgselect itemcget r$i -fill]
    if {$new != $old} {
      $w.imgselect itemconfigure t$i -text ""
      $w.imgselect itemconfigure r$i -fill $new
    }
    incr i
  }
}

proc ::ColorPicker::redrawpalette {} {
  variable w
  variable colors
  variable squaresize
  variable paletterows
  variable showlabels
  set hght [expr $squaresize*(ceil([llength $colors]/$paletterows)+1)]
  $w.frm2.imgpalette delete "all"
  $w.frm2.imgpalette configure -height [expr min($hght,500)] -scrollregion "0 0 0 $hght"
  set i 0
  set vis "hidden"
  if {$showlabels} {
    set vis "normal"
  }
  foreach col $colors {
    lassign $col nm r g b h s v
    set l [expr $i/$paletterows]
    set p [expr $i%$paletterows]
    $w.frm2.imgpalette create rectangle [expr $squaresize*$p] [expr $squaresize*$l] [expr $squaresize*($p+1)] [expr $squaresize*($l+1)] -fill [::ColorPicker::rgb2hex $r $g $b]
    set txtcol #000000
    if {$v<128} {set txtcol #ffffff}
    $w.frm2.imgpalette create text [expr int($squaresize*(0.5+$p))] [expr int($squaresize*(0.5+$l))] -text $nm -font "Helvetica, 6" -state $vis -tag t$i -fill $txtcol
    incr i
  }
}

proc ::ColorPicker::changepalette {} {
  set filename [tk_getOpenFile]
  if {[string length $filename]==0} {return}
  readcolors $filename
  redrawpalette
}

proc ::ColorPicker::addpalette {} {
  set filename [tk_getOpenFile]
  if {[string length $filename]==0} {return}
  readcolors $filename 1
  redrawpalette
}

proc ::ColorPicker::sortpalette {} {
  variable colors
  set colors [lsort -integer -index 4 [lsort -integer -index 6 [lsort -integer -index 5 $colors]]]
  redrawpalette
}

proc ::ColorPicker::selectcolor {x y} {
  variable w
  variable selectorrows
  variable rectwdth
  variable recthght
  variable texthght
  variable selected
  set x [expr int([$w.imgselect canvasx $x]/$rectwdth)]
  set y [expr int([$w.imgselect canvasy $y]/($texthght+$recthght))]
  set i [expr $selectorrows*$y+$x]
  if {$i >= [colorinfo num]} { return }
  set x [expr $rectwdth*($i%$selectorrows)]
  set y [expr ($texthght+$recthght)*($i/$selectorrows)]
  set p [expr $i%$selectorrows]
  set l [expr $i/$selectorrows]
  set selected $i
  $w.imgselect coords marker [expr $rectwdth*$p] [expr ($texthght+$recthght)*$l] [expr $rectwdth*($p+1)] [expr ($texthght+$recthght)*($l+1)]
#   after 500
  update
}

proc ::ColorPicker::setcolor {x y} {
  variable w
  variable paletterows
  variable squaresize
  variable colors
  variable colorsort
  variable selected
  set x [expr int([$w.frm2.imgpalette canvasx $x]/$squaresize)]
  set y [expr int([$w.frm2.imgpalette canvasy $y]/$squaresize)]
  set i [expr $paletterows*$y+$x]
  if {$i >= [llength $colors]} { return }
  lassign [lindex $colors $i] nm r g b h s v
  set r [expr $r/255.0]
  set g [expr $g/255.0]
  set b [expr $b/255.0]
  set txtcol #000000
  if {$v<128} {set txtcol #ffffff}
  color change rgb [lindex $colorsort $selected] $r $g $b
  $w.imgselect itemconfigure r$selected -fill [vmd2hex $r $g $b]
  $w.imgselect itemconfigure t$selected -fill $txtcol -text $nm
#   after 500
  update
}

proc ::ColorPicker::printcolorscript {} {
  set filename [tk_getSaveFile]
  if {[string length $filename]==0} {return}
  set fp [open $filename w]
  for {set i 0} {$i<[colorinfo num]} {incr i} {
    lassign [colorinfo rgb $i] r g b
    puts $fp [format "color change rgb %2d %8.6f %8.6f %8.6f" $i $r $g $b]
  }
  close $fp
}

proc ::ColorPicker::setlabels {} {
  variable showlabels
  variable colors
  variable w
  set vis "hidden"
  if {$showlabels} { set vis "normal" }
  set mx [llength $colors]
  for {set i 0} {$i < $mx} {incr i} {
      $w.frm2.imgpalette itemconfigure t$i -state $vis
  }
  for {set i 0} {$i<[colorinfo num]} {incr i} {
    $w.imgselect itemconfigure t$i -state $vis
  }
}

proc colorpicker_tk_cb {} {
  ::ColorPicker::colorpicker
  return $::ColorPicker::w
}