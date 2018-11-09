#!/usr/bin/perl

@files = <mpi*.e*>;

$nfiles = @files;

for($i=0;$i<$nfiles;$i++) {
    $files[$i] =~ /mpitest_([1234567890]+)_([1234567890]+)\.e/;
    $n = $1;
    $p = $2;
    open FILE, "<$files[$i]" or die$!;
    while ($line = <FILE>) {
        chomp $line;
        if ($line =~ /^real[ \t]+(.*)/) {
            $line = $1;
            $line =~ /([0123456789\.]+)m([0123456789\.]+)s/;
            $time = $1*60.0+$2;
            print "$n\t$p\t$time\n";
        }
    }
    close FILE;
}


