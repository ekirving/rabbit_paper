use File::Basename;
my $files = $ARGV[0];
my $fastq = $ARGV[1];
open MY, "$files";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

my %allBam;
my %countBam;
while (my $line = <MY>) {
	chomp($line);
	my @w = split (' ', $line); 
	
	#GET ID + LIB
	my $lib = $w[1];
	my $id = $w[0];
#	system ("rm $id/*.gz"); #WARNING COMMENT OUT
#	system ("rm $id/*.bam");
	system ("mkdir $id");

	my $shell = $id.'_'.$mon.'_'.$year.'.sh';
	if (exists $allBam{$id}) {
		open SH, ">>$shell";
	}	
		
	else {
		open SH, ">$shell";
	}
	
	#GET PATH FASTQ 
#	($name,$dir,$ext) = fileparse($line,'\..*');
	open IN, "cat $fastq | grep $lib | grep $id | grep fastq |";
	my $count=0;
	my %fastq;
	
	while (my $l = <IN>) {
		++$count;
		chomp($l);
		$fastq{$count} = $l;
	}
	
	my $fastq_toAln;
	foreach my $c (keys %fastq) {
		#ADRM 
		#/usr/local/bin/AdapterRemoval --file1 $i --trimns --output1 `basename $i .fq.gz`_adRm.fq.gz --gzip --basename `basename $i .fq.gz		
		($name,$dir,$ext) = fileparse($fastq{$c},'\..*');
		my @s = split ('_', $name);
		my $f = $s[$#s];
		my $l = $#s-2;
		my $lane=$s[$l];
		#TOG_R3YV_KD033_27_TATGCA_L003_R1_066.fastq.gz
		my $adRm_out = "$id".'/'."$id".'_'.$lib.'_'.$lane.'_'.$f;
		my $tmp_f =  $adRm_out = "$id".'/'."$id".'_'.$lib.'_'.$lane.'_'.$f."-adRm.fq.gz";
		my $cmd_adRm_out = "/usr/local/bin/AdapterRemoval --file1 $fastq{$c} --trimns --output1 $tmp_f --gzip --basename $adRm_out";
		unless (-e $tmp_f) {
			print SH "$cmd_adRm_out\n"; 
			$fastq_toAln.=" $tmp_f";
		}
	
	}
	
	#Do alignment from STDIN
	my $sai = "$id".'/'."$id".'_'.$lib.'.sai';
	my $bam = "$id".'/'."$id".'_'.$lib.'.bam';  	
	my $log = "$id".'/'."$id".'_'.$lib.'.log';
	unless (-e $bam) {
		my $cmd_aln = "gunzip -c $fastq_toAln | bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa - 1> $sai 2> $log";
		my $cmd_aln1 = "gunzip -c $fastq_toAln | bwa samse -r \'\@RG\\tID:$lib\\tLB:$lib\\tPL:ILLUMINA\\tSM:$id\' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa $sai - | samtools view -Shu - 1> $bam 2>> $log";
		print SH "$cmd_aln\n$cmd_aln1\n";
	}
	$allBam{$id} .= " $bam";
	++$countBam{$id};
}

foreach my $id (keys %allBam) {
	my $sort = "$id".'/'."$id".'_sorted';
	my $sorted = "$id".'/'."$id".'_sorted.bam';
	my $shell = $id.'_'.$mon.'_'.$year.'.sh';
        open SH, ">>$shell";
	if ($countBam{$id} > 1) {
		my $bam_merged="$id".'/'."$id".'_merged.bam';
		print SH "samtools merge - $allBam{$id} > $bam_merged\n";
		print SH "samtools index $bam_merged\n";
		print SH "samtools sort $bam_merged $sort\nsamtools index $sorted\n";	
	}
	else {
		 print SH "samtools sort $allBam{$id} $sort\nsamtools index $sorted\n";			
	}

	$RC="$id".'/'."$id".'_RC.txt';

	$rmdup = "$id".'/'."$id".'_rmdup.bam';
	$dup = "$id".'/'."$id".'_dup.txt';
	print SH "samtools rmdup -s $sorted $rmdup 2>&1 >/dev/null | cut -d' ' -f 6 >  $dup\n";
	print SH "samtools index $rmdup\n";

	print SH "samtools view $rmdup | perl /media/raid5/laurent/full_run_results/bin/compute_l.pl > $RC\n";
	
	$cov = "$id".'/'."$id".'.cov';
	print SH "samtools view -q 20 -b $rmdup | bedtools genomecov -ibam stdin | grep -v genome > $cov\n";
	
	
}
