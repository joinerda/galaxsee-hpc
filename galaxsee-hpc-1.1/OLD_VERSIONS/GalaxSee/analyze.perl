#!/usr/bin/perl

$j=0;
while($line=<STDIN>) {
    if($line=~/n = (.*)/) {
        $n=$1;
    }
    if($line=~/(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)/) {
        $i = $2;
        if($i==1) {
            $time[$j]=$1;
        }
        $x[$j][$i]=$3;
        $y[$j][$i]=$4;
        $z[$j][$i]=$5;
        $vx[$j][$i]=$6;
        $vy[$j][$i]=$7;
        $vz[$j][$i]=$8;
        if($i==$n-1) {
            $j++;
        }
    }
}
$ntime = $j;
print "n = $n\tntime = $ntime\n";
for($j=0;$j<$ntime;$j++) {
    $rave=0.0;
    $vave=0.0;
    for($i=0;$i<$n;$i++) {
        $r = sqrt($x[$j][$i]*$x[$j][$i]+
                  $y[$j][$i]*$y[$j][$i]+
                  $z[$j][$i]*$z[$j][$i]);
        $rave = $rave + $r;
        $v = sqrt($vx[$j][$i]*$vx[$j][$i]+
                  $vy[$j][$i]*$vy[$j][$i]+
                  $vz[$j][$i]*$vz[$j][$i]);
        $vave = $vave + $v;
    }
    $rave = $rave/$n;
    $vave = $vave/$n;
    print "$time[$j]\t$rave\t$vave\n";
}
