#!/usr/bin/perl -w
use strict;

# Compare2.pl
#
# A script to compare the ordering results produced by my Hi-C-based method (Lachesis) and Adam Waalkes' RAD tag-based method.
#
# Shares some code with Compare2.pl, which compares clustering results.
#
# Josh Burton
# August 2013



# INPUT PARAMETERS
my $trunk_only = 0;
my $min_Q = 0; # if this goes above 21, the tags don't match up... grr
my $include_scaffold_name = 1; # this option adds the scaffold name to the dotplots so they're more human-readable
my $include_scaffold_length = 1; # this option adds the scaffold length to the dotplots for extra plotting goodness
#my $case = "20"; # "20", "10", "20_split" - TEMP - master control
my $case = "10"; # TEMP


# INPUT FILES
my $assembly_fasta = "$ENV{'HOME'}/L/rhododendron/fragScaff/final.assembly.fasta";
my $scaffold_names_file    = "$assembly_fasta.names";
my $scaffold_lengths_file  = "$assembly_fasta.dict";
#my $scaffold_lengths_file2 = "$ENV{'HOME'}/L/rhododendron/assembly_nosplit/assembly.fasta.dict";
my $RAD_tags_clusters_file = "from_Adam/latest.txt";
#$RAD_tags_clusters_file = "from_Adam/final_order_split.txt" if $case =~ /split/;
my $out_dir = "$ENV{'HOME'}/L/git/out/rhododendron_mod";
#$out_dir .= "_1010" if $case =~ /10/;
#$out_dir .= "_split" if $case =~ /split/;
my @Lachesis_ordering_files = $trunk_only ?
    glob "$out_dir/cached_data/group*.trunk.ordering" :
    glob "$out_dir/main_results/group*.ordering";

my %scaffold_lengths; # will be filled by &get_scaffold_lengths

my @LG_names = ( "LG01", "LG02", "LG03", "LG04", "LG05", "LG06", "LG07", "LG08", "LG09", "LG10", "LG11", "LG12", "LG13" );

# Set of linkage groups to manually flip.
my %LGs_to_reverse = ( 'LG01' => 1, 'LG02' => 1, 'LG03' => 1, 'LG06' => 1, 'LG10' => 1, 'LG12' => 1 );




# Copied from http://www.perlmonks.org/?node_id=435711
sub commify {
    my ( $sign, $int, $frac ) = ( $_[0] =~ /^([+-]?)(\d*)(.*)/ );
    my $commified = (
        scalar reverse join ',',
        unpack '(A3)*',
        scalar reverse $int
    );
    return $sign . $commified . $frac;
}





# Sorting function that takes subscript letters into account.  E.g., 11 < 12 < 12a < 12b < 13.
sub scaffold_sort {
    $a =~ /(\d+)(\w*?)/;
    my ($a1, $a2) = ($1,$2);
    $b =~ /(\d+)(\w*?)/;
    my ($b1, $b2) = ($1,$2);
    #print "STUFF: $a,$b,$a1,$b1,$a2,$b2\n";
    return $a1 <=> $b1 || $a2 cmp $b2;
}


# get_scaffold_names: Read scaffold names from a *.fasta.names file.  Return an array of scaffold names, and a hash of the scaffolds that have been split.
sub get_scaffold_names($) {
    
    open IN, '<', $_[0] or die "Can't find scaffold fasta names file $_[0]: $!\n";
    
    my @names = ();
    my %splits = ();
    
    
    foreach (<IN>) {
	
	die unless /^scaffold_([\w\.]+)\r?$/;
	my $name = $1;
	push @names, $name;
	
	# Additionally, record which scaffolds have been split (e.g., "4b".)
	if ( $name =~ /(\d+)(\D+)/ ) {
	    $splits{$name} = 1;
	    $splits{$1} = 1;
	}
    }
    
    close IN;
    
    return (\@names, \%splits );
}


my %tag_to_LG;


# get_R_orderings: Get the groups out of a text file.  Return a hash of <lowest contig ID> -> <array of contigs, in order>.
sub get_R_orderings($ )
{
    my %orderings = (); # this will contain output
    
    my $verbose = 1;
    
    my ($file) = @_;
    die "Can't find file: $file\n" unless -e $file;
    print "get_R_orderings from file: $file\n" if $verbose;
    
    # Loop over all lines in the file.
    open FILE, '<', $file or die;
    
    my $LG_counter = 0;
    
    foreach (<FILE>) {
	
	next if /^\#/; # skip commented lines
	
	# Name this linkage group.  This name goes on the output files.
	my $LG = $LG_names[$LG_counter++];
	
	# Get the contig names.
	my @tokens = split;
	
	# Find the lowest element of the tokens that is NOT a split scaffold.  Split scaffolds are unreliable because they're placed twice in the RAD tags.
	my @tokens_sorted = sort scaffold_sort @tokens;
	my $i = 0;
	while ( $tokens_sorted[$i] =~ /[^0-9]/ ) { $i++ } # split scaffolds contain non-numeric characeters e.g., "0.xx", "4a"
	
	
	# Use this scaffold ID as a tag to describe this cluster.
	my $tag = $tokens_sorted[$i];
	
	#if ( $trunk_only ) { # deal with the fact that not all contigs appear in the trunks
	#    $tag = 70  if $tag == 0;
	#    $tag = 17  if $tag == 1;
	#    $tag = 9   if $tag == 2;
	#    $tag = 20  if $tag == 3;
	#    $tag = 15  if $tag == 8;
	#    $tag = 62  if $tag == 13;
	#    $tag = 149 if $tag == 23;
	#}
	
	$tag_to_LG{$tag} = $LG;
	print "LG = $LG\tTAG = $tag\n" if $verbose;
	
	$orderings{$tag} = \@tokens;
    }
    
    close FILE;
    print "\n" if $verbose;
    
    
    return \%orderings;
}


# get_L_orderings: Get the orderings out of a set of text files, i.e., group*.ordering, produced by Lachesis.
# Return a hash of <lowest contig ID> -> <array of contigs, in order>.
sub get_L_orderings($@ )
{
    my ($splits, @files) = @_;
    
    my %orderings = (); # this will contain output
    
    my $verbose = 0;
    
    # Loop over all files.  For each file, compile an array containing the ordering.
    foreach my $file (@files) {
	die "Can't find file: $file\n" unless -e $file;
	print "get_L_orderings from file: $file\n";
	
	my @ordering = ();
	
	# Loop over all lines in the file.
	open FILE, '<', $file or die;
	
	foreach (<FILE>) {
	    
	    next if /^\#/; # skip commented lines
	    
	    # The second token in the file should be a contig name.  Get the ID out of the name.
	    my @tokens = split;
	    die unless $tokens[1] =~ /scaffold_([\w\.]+)/;
	    next unless $tokens[3] >= $min_Q;
	    push @ordering, $1;
	}
	
	close FILE;
	next unless @ordering; # skip orderings with no data (e.g., singleton clusters)
	
	# Find the lowest element of the tokens that is NOT a split scaffold.  Split scaffolds are unreliable because they're placed twice in the RAD tags.
	my @ordering_sorted = sort scaffold_sort @ordering;
	my $i = 0;
	while ( $splits->{ $ordering_sorted[$i] } ) { $i++ }
	
	# Use this scaffold ID as a tag to describe this cluster.
	my $tag = $ordering_sorted[$i];
	print "$tag\t->\t@ordering_sorted\n" if $verbose;
	
	$orderings{$tag} = \@ordering;
    }
   
    return \%orderings;
}



# get_scaffold_lengths: Open the input filenames and fill the hash %scaffold_lengths.
sub get_scaffold_lengths(@) {
    
    %scaffold_lengths = ();
    
    foreach my $infile (@_) {
    
    # Read in both scaffold lengths files, to get both the split and non-split versions of split scaffolds.
	open IN, '<', $infile or die "get_scaffold_lengths: Can't find input file $infile: $!";
	foreach (<IN>) {
	    /^\@SQ\s+SN\:scaffold_([\w\.]+)\s+LN\:(\d+)\s+/ or next;
	    $scaffold_lengths{$1} = $2;
	}
	close IN;
    }
}



# report_on_scaffolds: Report the number and total length of scaffolds listed in an array.
sub report_on_scaffolds($@) {
    
    my ($label, @scaffolds) = @_;
    my $N_scaffolds = scalar @scaffolds;
    
    my $total_len = 0;
    foreach my $scaffold (@scaffolds) {
	$total_len += $scaffold_lengths{$scaffold};
    }
    $total_len = commify $total_len;
    
    print "$N_scaffolds scaffolds, with total length $total_len: $label\n";
}





print "STARTING: \$case = $case, \$trunk_only = $trunk_only, \$min_Q = $min_Q\n";
print "Number of Lachesis orderings: ", scalar @Lachesis_ordering_files, "\n";

my ($scaffold_names, $split_scaffolds) = get_scaffold_names( $scaffold_names_file );

# Read in both scaffold lengths files, to get both the split and non-split versions of split scaffolds.
&get_scaffold_lengths($scaffold_lengths_file);


my $RAD_tag_orderings  = &get_R_orderings( $RAD_tags_clusters_file );
my $Lachesis_orderings = &get_L_orderings( $split_scaffolds, @Lachesis_ordering_files );

# Verify that the two sets of clustering results have the same set of tags.
# If this assert fails, the groups in the two sets of clustering results probably don't match up nicely.  (Have I been tweaking the Lachesis results?)
my @tags  = sort scaffold_sort keys %$Lachesis_orderings;
my @tags2 = sort scaffold_sort keys %$RAD_tag_orderings;
print "Lachesis group tags: (", scalar @tags , ")\t@tags\n";
print "RAD-tags group tags: (", scalar @tags2, ")\t@tags2\n";
map { die unless $tags[$_] eq $tags2[$_] } ( 0..$#tags );



# Unresolved ambiguities in Adam's RAD-tag linkage map.
# If a set of N consecutive scaffolds can't be placed in order w.r.t. one another, *all but the last of them* should be listed here.
# The values indicate the linkage group the scaffolds belong in, and aren't actually used.
my %RAD_tag_ambiguities = ();
#my %RAD_tag_ambiguities = ( 458 => 1, 4312 => 1, 5329 => 1, 4307 => 2, 1391 => 2, 513 => 2, 526 => 2, 1822 => 2, 914 => 2, 1023 => 2, 1089 => 2, 675 => 2, 1175 => 2, 1429 => 2, 1765 => 2, 2625 => 2, 2889 => 2, 3159 => 2, 4582 => 2, 4271 => 2, 1328 => 2, 1479 => 2, 3152 => 2, 1588 => 3, 1393 => 3, 1447 => 3, 1542 => 3, 1366 => 3, 1553 => 3, 3366 => 3, 433 => 3, 838 => 3, 2265 => 3, 2344 => 3, 961 => 3, 1151 => 4, 922 => 4, 1653 => 4, 624 => 4, 2115 => 5, 934 => 5, 6575 => 5, 1297 => 5, 294 => 5, 928 => 5, 454 => 5, 1885 => 5, 1798 => 5, 1702 => 5, 2126 => 5, 2317 => 6, 2087 => 6, 2269 => 6, 233 => 6, 559 => 6, 580 => 6, 620 => 6, 1153 => 6, 1418 => 6, 5041 => 6, 886 => 6, 1267 => 6, 2580 => 6, 645 => 6, 1338 => 6, 1485 => 6, 1916 => 6, 1209 => 6, 2290 => 6, 3590 => 6, 3972 => 6, 2080 => 6, 1886 => 6, 1941 => 6, 1217 => 7, 1360 => 7, 2689 => 7, 2182 => 7, 717 => 7, 888 => 7, 1948 => 7, 796 => 7, 899 => 7, 1246 => 7, 1710 => 7, 2305 => 7, 2757 => 7, 865 => 7, 1905 => 7, 977 => 8, 2918 => 8, 4766 => 8, 2641 => 8, 2253 => 9, 3091 => 9, 1042 => 9, 6539 => 9, 1643 => 9, 1252 => 9, 1429 => 9, 1908 => 9, 2177 => 9, 3216 => 9, 857 => 9, 2102 => 9, 405 => 9, 564 => 9, 1325 => 9, 973 => 9, 1280 => 9, 1446 => 9, 1013 => 10, 1752 => 10, 1143 => 10, 1501 => 10, 2861 => 10, 1309 => 10, 1475 => 10, 1845 => 10, 1917 => 10, 3118 => 10, 1947 => 10, 984 => 11, 2202 => 11, 499 => 11, 852 => 11, 1099 => 11, 731 => 11, 754 => 11, 1072 => 11, 2006 => 11, 5312 => 11, 900 => 11, 1294 => 12, 283 => 12, 1040 => 12, 1434 => 12, 291 => 12, 2327 => 12, 3146 => 12, 4331 => 12, 2284 => 12, 2082 => 12, 1540 => 12, 1833 => 12, 2093 => 12, 1678 => 13, 1354 => 13, 2191 => 13, 632 => 13, 943 => 13, 1416 => 13, 1537 => 13, 1754 => 13, 1759 => 13, 1279 => 13, 1743 => 13, 1120 => 13, 1484 => 13, 1028 => 13, 1209 => 13, 1521 => 13, 1086 => 13, 1490 => 13 );



my ( @in_both, @L_only, %R_only );

# Now, make some awesome dotplots.
# For each tag (i.e., each group), find contigs that are placed in both orderings.  Each contig's positions in the two orderings imply an ordered pair (x,y).
# The dotplot is simply a plot of all such ordered pairs.
foreach my $tag (@tags) {
    #next unless $tag == 1; # TEMP: group1 (Adam's LGO1)
    
    my $file_label = "LR_comparison.20.20";
    $file_label = "LR_comparison.10.10" if $case eq '10';
    $file_label = "LR_comparison.20.20.split" if $case eq '20_split';
    my $outfile = "$file_label.$tag_to_LG{$tag}";
    $outfile .= ".trunk" if $trunk_only;
    $outfile .= ".Q$min_Q" if $min_Q;
    print "Writing to file $outfile...";
    open OUT,  '>', $outfile or die;
    open OUT2, '>', $outfile . ".2" or die;
    print OUT2 "chrom\tx\txend\ty\tyend\torient_error\n";
    
    # Make lookup tables of [contig name] -> [position in the Lachesis and RAD-tag orderings].
    my (%contig_to_Lpos, %contig_to_Rpos);
    
    my @L_ordering = @{$Lachesis_orderings->{$tag}};
    my @R_ordering = @{$RAD_tag_orderings ->{$tag}};
    #print "\n@L_ordering\n@R_ordering\n";
    my $i = 0;
    foreach (@L_ordering) {
	$contig_to_Lpos{$_} = $i;
	$i++;
    }
    $i = 0;
    foreach (@R_ordering) {
	$contig_to_Rpos{$_} = $i;
	$i++ unless $RAD_tag_ambiguities{$_}; # don't increment certain RAD-tag orderings - corresponding to unresolved ambiguities in Adam's linkage map
    }
    
    # Find all contigs that appear in both orderings, and find their ordered pairs (x,y), where x = position in RAD tags, y = position in Lachesis.
    # Also for each contig, print two more numbers: the scaffold name (which affects color plotting in QuickDotplot) and length (also useful for plotting!)
    my %scaff_points = ();
    foreach my $contig ( sort scaffold_sort keys %contig_to_Lpos ) {
	my $contig2 = $contig;
	$contig2 =~ s/\D$//; # remove the subscript "a/b" from Lachesis' contig names for direct comparison with RAD tag contig names
	unless ( exists $contig_to_Rpos{$contig2} ) { # contig must appear in both orderings to be plotted
	    push @L_only, $contig;
	    next;
	}
	push @in_both, $contig2;
	
	my $x = $contig_to_Rpos{$contig2};
	my $y = $contig_to_Lpos{$contig};
	my @point = ($x,$y);
	push @point, "scaffold_$contig" if $include_scaffold_name;
	push @point, $scaffold_lengths{$contig} if $include_scaffold_length;
	$scaff_points{$contig} = \@point;
    }
    
    # Bookkeeping.
    map { $R_only{$_}++ unless ( exists $contig_to_Lpos{$_} ) } keys %contig_to_Rpos;
    
    
    # For each of the orderings (Lachesis, RAD-tag), get the contig ordering.  Flip to reorient if necessary.
    my @contigs_in_L_order = sort {$contig_to_Lpos{$a} <=> $contig_to_Lpos{$b} } keys %contig_to_Lpos;
    my @contigs_in_R_order = sort {$contig_to_Rpos{$a} <=> $contig_to_Rpos{$b} } keys %contig_to_Rpos;
    @contigs_in_R_order = reverse @contigs_in_R_order if $LGs_to_reverse{$tag_to_LG{$tag}};
    
    # Now, for each ordering, make lookup tables of [position in ordering] -> [position *BY BP* in ordering].
    my ($Lseq, $Rseq) = (0,0);
    my (%contig_to_Lseq, %contig_to_Rseq);
    foreach my $contig ( @contigs_in_L_order ) {
	$contig_to_Lseq{$contig} = $Lseq;
	$Lseq += $scaffold_lengths{$contig};
    }
    foreach my $contig ( @contigs_in_R_order ) {
	$contig_to_Rseq{$contig} = $Rseq;
	$Rseq += $scaffold_lengths{$contig};
    }
    
    
    
    foreach my $contig ( @contigs_in_L_order ) { # TEMP: or just sort scaffold_sort keys %contig_to_Lpos
	print OUT join ("\t", @{$scaff_points{$contig}} ), "\n"
	    if exists $scaff_points{$contig};
    }
    
    
    # Print the file that shows the contigs as slanty diagonal lines instead of dots.  This is for use with Lineplot.R.
    # NOTE: For rhododendron, x_start and x_stop don't make much sense because there's no orientation info in Adam's LGs.  The Lineplot plotting reflects this.
    foreach my $contig (@in_both) {
	next unless exists $contig_to_Lseq{$contig}; # TEMP: not sure why this is necessary, it shouldn't be, but I don't care
	my $x_start = $contig_to_Rseq{$contig};
	my $y_start = $contig_to_Lseq{$contig};
	my $x_stop  = $x_start + $scaffold_lengths{$contig};
	my $y_stop  = $y_start + $scaffold_lengths{$contig};
	print OUT2 "chrom\t$x_start\t$x_stop\t$y_start\t$y_stop\t0\n";
    }
    
    
    close OUT;
    close OUT2;
    
    print "\t", scalar keys %scaff_points, " points written.\n";
    
    # Now use QuickDotplot to convert the dotplot to a jpeg image.
    my $legend_flag = $include_scaffold_name ? "-legend" : "";
    system ( "QuickDotplot $legend_flag $outfile" );
    #system ( "rm $outfile" );
}


@L_only = sort scaffold_sort @L_only;
my @R_only = ();
map { push @R_only, $_ if $R_only{$_} == 1 } sort scaffold_sort keys %R_only;

report_on_scaffolds( "in both sets of orderings", @in_both );
report_on_scaffolds( "in Lachesis orderings only", @L_only );
report_on_scaffolds( "in RAD-tag orderings only", @R_only );
