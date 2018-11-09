#!/usr/bin/perl

@n = qw(500 1000 2000 5000);
@procs = qw(1 2 3 4 8 16 32 64);

$n_n = @n;
$n_procs = @procs;

for($i=0;$i<$n_n;$i++) {
    for($j=0;$j<$n_procs;$j++) {
        $file = "mpitest_$n[$i]_$procs[$j].qsub";
        open(QSUB,"> $file");
        print QSUB "#PBS -N mpitest_$n[$i]_$procs[$j]\n";
        print QSUB "#PBS -q default\n";
        $num_nodes = int($procs[$j]/8);
        $num_cores = $procs[$j]%8;
        if($num_cores != 0) {$num_nodes=1;}
        if($num_cores == 0) {$num_cores=8;}
        $num_procs = $num_nodes*$num_cores;
        print QSUB "#PBS -l nodes=$num_nodes:ppn=$num_cores\n";
        print QSUB "#PBS -l cput=10:00:00\n";
        print QSUB "PROGRAM='./galaxsee n$n[$i].gal'\n";
        print QSUB "cd \$PBS_O_WORKDIR\n";
        print QSUB "time /opt/mpich/intel/bin/mpirun -machinefile \$PBS_NODEFILE -np $num_procs \$PROGRAM\n";
        print QSUB "exit \n";
        `qsub $file`;
        close(QSUB);
    }
}


