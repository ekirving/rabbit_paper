my $rep = $ARGV[0];
my $prefix = $ARGV[1];
my $dir = $ARGV[2];
my %column;
my $count;
my %hash;
my $ncol;
my %head;
while (my $line = <STDIN>) {
        chomp($line);
        ++$count;
        my @w = split ("\t", $line);
        if ($count==1) {
                my $col = $#w - 6;
                $ncol = $#w;
                for my $u (1 .. $rep) {
                        for my $z (0 .. $col) {
                                my $rand = int(rand($col));
                                my $colX = $rand+6;
                                my $cc = $z + 6;
                                $column{$u}{$cc} = $colX;
                        }
                }
        }
        for my $i (0 .. 5) {
                $head{$w[1]} .= $w[$i]."\t";
        }
        for my $i (6 .. $#w) {
                $hash{$w[1]}{$i} = "$w[$i]";
        }
}
for my $i (1 .. $rep) {
        my $ped = $dir.$prefix."_".$i.'.ped';
        my $map = $dir.$prefix."_".$i.'.map';
        my $mapo = $prefix.'.map';
        system ("cp $mapo $map");
        open OUT, ">$ped";
        foreach my $key (keys %head) {
                print OUT "$head{$key}";
                for my $col (sort {$a <=> $b} keys %{ $hash{$key} }) {
                        print OUT "$hash{$key}{$column{$i}{$col}}\t";
                }
                print OUT "\n";
        }
}