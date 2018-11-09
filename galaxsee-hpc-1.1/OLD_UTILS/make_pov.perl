#!/usr/bin/perl

$infile = $ARGV[0];
if ($infile =~ /(.*)\.dump$/) {
    open (FILE,"<$infile") || die "Cannot open file";
    $base = $1;
    $outfile = $base.".pov";
    open (OUT,">$outfile") || die "Cannot open file";
} else {
    die "no file given";
}

$angle = $base;
$angle =~ /([1234567890]+)/;
$angle = $1/(90101);
$angle = $angle*2.0*3.14159;

print $angle."\n";
$xpos = 4.0*sin($angle+1.0);
$zpos = 4.0*cos($angle+1.0);
 

print OUT "#include \"colors.inc\"\n";
print OUT "#include \"myShapes.inc\"\n";
print OUT "#include \"basicView.inc\"\n";
print OUT "basicView(<$xpos,1.7,$zpos>,<0,0,0>,45)\n";

print OUT "#local Color = color rgb<1.0,1.0,1.0>;\n";
$line = <FILE>;
    chomp $line;
($n,$time,$scale,$rest) = split(/[\t ]/,$line);
#print OUT "text { ttf \"timrom.ttf\" \"t = $time\" 1, 0 pigment { White } }\n";
while($line=<FILE>) {
    chomp $line;
    ($i,$x,$y,$z,$rest) = split(/[\t ]/,$line);
    $x = $x/$scale;
    $y = $y/$scale;
    $z = $z/$scale;
    print OUT "point(<$x,$y,$z>,0.1,0.9,Color)\n";
    print OUT "point(<$x,$y,$z>,0.05,0.6,Color)\n";
    print OUT "point(<$x,$y,$z>,0.02,0.2,Color)\n";
}
close(FILE);
