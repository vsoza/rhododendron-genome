#!/usr/bin/perl -w
use strict;

# Compare.pl
#
# A script to compare the clustering results produced by my Hi-C-based method (Lachesis) and Adam Waalkes' RAD tag-based method.
# Also writes a merged clusters file (clusters.merged_with_RAD.by_name.txt) with the permissive union of the two clustering results.
#
# Like Compare.py, but I'm tired of Python.  Shares some code with Compare2.pl, which compares ordering results.
#
# Josh Burton
# August 2013


# INPUT FILES
my $scaffold_names_file   = "$ENV{'HOME'}/L/rhododendron/fragScaff/final.assembly.fasta.names";
my $scaffold_lengths_file = "$ENV{'HOME'}/L/rhododendron/fragScaff/final.assembly.fasta.FastaSize";
my $RAD_tags_clusters_file = "from_Adam/latest.txt";
my $Lachesis_dir = "$ENV{'HOME'}/L/git/out/rhododendron";
my $Lachesis_clusters_file = "$Lachesis_dir/main_results/clusters.by_name.txt";
# OUTPUT FILE
my   $merged_clusters_file = "$Lachesis_dir/main_results/clusters.merged_with_RAD.txt";



# Count the number of lines in the names file to get the number of contigs.
open IN, '<', $scaffold_names_file; while (<IN>) {};
my $N_contigs = $.;


my @comments;



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
    
    
    while (<IN>) {
	
	# Older form of scaffold names: preceded by "original_scaffold_" or "fragScaff_scaffold_"
	#die unless /^(original|fragScaff)_scaffold_(\w+)\r?$/;
	#my $name = $2;
	die "wtf: $_" unless /^(scaffold_)?([\w\.]+)\r?$/;
	my $name = $2;
	#$name =~ s/scaffold_//;
	push @names, $name;
    }
    
    close IN;
    
    die "Inconsistent value of \$N_contigs: $N_contigs vs. ", scalar @names unless scalar @names == $N_contigs; # if this fails, change hard-wired $N_contigs
    
    return \@names;
}



# by_scaffold_ID: Helper function for get_clusters' scaffold sorting.
sub by_scaffold_ID {
    my $ID1 = $a;
    my $ID2 = $b;
    
    # Trim scaffold IDs off.
    $ID1 =~ s/scaffold_//;
    $ID2 =~ s/scaffold_//;
    
    # Account for the possible suffixes, i.e., 4a, 4b.
    my $ID1_int = $ID1;
    my $ID2_int = $ID2;
    $ID1_int =~ s/\D$//;
    $ID2_int =~ s/\D$//;
    
    return $ID1_int <=> $ID2_int || $ID1 cmp $ID2;
}



# TEMP: list of scaffolds that don't match between the Lachesis and RAD tag assemblies, so don't use them as tags
my @is_bad = (1);
my %is_bad; # convert array to hash
map { $is_bad{$_} = 1 } @is_bad;


# get_clusters: Get the groups out of a text file.  Return a hash of <lowest unique contig ID> -> <array of contigs, in order>.
# $is_RAD_tags: Boolean flag indicating whether these are Adam's RAD tag results instead of mine.
sub get_clusters($$)
{
    my %clusters = (); # this will contain output
    
    my $verbose = 0;
    
    my ($file,$is_RAD_tags) = @_;
    die "Can't find file: $file\n" unless -e $file;
    print "get_clusters from file: $file\n" if $verbose;
    
    
    
    
    # Now, for reals.
    # Loop over all lines in the file.
    open FILE, '<', $file or die;
    
    while (<FILE>) {
	
	# Skip commented lines, but record them for use later.
	if ( /^\#/ ) { push @comments, $_; next; }
	
	# Get the contig names.  For Lachesis results, parse out the "scaffold_".
	my @tokens_unsorted = split;
	next if scalar @tokens_unsorted < 2; # skip singleton clusters and blank lines
	map { s/scaffold_// } @tokens_unsorted;
	
	# Find the lowest element of the tokens that is NOT a split scaffold.  Split scaffolds are unreliable because they're placed twice in the RAD tags.
	my @tokens = sort by_scaffold_ID @tokens_unsorted;
	my $i = 0;
	while ( $tokens[$i] =~ /\D/ )  { $i++ } # don't use weird scaffold names, i.e. "0.x" "4b" as representatives of their clusters
	while ( $is_bad{$tokens[$i]} ) { $i++ } # these are bad
	#unless ($is_good{$tokens[$i]} ) { $i++ } # TEMP: these are the tokens I've chosen manually to use
	
	# Use this scaffold ID as a tag to describe this cluster.
	my $tag = $tokens[$i];
	#$tag = 32 if $tag == 19; # TEMP: covering up for a weird bug that I don't wanna fix
	#print "$is_RAD_tags, $i -> $tag\n"; 
	
	print "$tag\t->\t@tokens\n" if $verbose;
	
	$clusters{$tag} = \@tokens_unsorted;
    }
    
    close FILE;
    print "\n" if $verbose;
    
    
    return \%clusters;
}


# report_on_scaffolds: Report the number and total length of scaffolds listed in an array.
sub report_on_scaffolds($@) {
    
    my ($label, @scaffolds) = @_;
    my $N_scaffolds = scalar @scaffolds;
    
    my %scaffold_lens;
    open IN, '<', $scaffold_lengths_file or die "Can't open scaffold lengths file $scaffold_lengths_file";
    while (<IN>) {
	# old fasta.dict format regex
	#/^\@SQ\s+SN\:scaffold_(\w+)\s+LN\:(\d+)\s+/ or next;
	#$scaffold_lens{$1} = $2;
	/^\s+(\d+).+scaffold_([\w\.]+)/ or next; # this filters out the 'TOTAL' line at the end
	$scaffold_lens{$2} = $1;
    }
    close IN;
    
    my $total_len = 0;
    foreach my $scaffold (@scaffolds) {
	$total_len += $scaffold_lens{$scaffold};
    }
    $total_len = commify $total_len;
    
    print "$N_scaffolds scaffolds, with total length $total_len: $label\n";
}







my ($scaffold_names) = &get_scaffold_names( $scaffold_names_file );


my $Lachesis_clusters = &get_clusters( $Lachesis_clusters_file, 0 );
my $RAD_tag_clusters  = &get_clusters( $RAD_tags_clusters_file, 1 );

# Verify that the two sets of clustering results have the same set of tags.
# If this assert fails, the groups in the two sets of clustering results probably don't match up nicely.  (Have I been tweaking the Lachesis results?)
my @tags  = sort scaffold_sort keys %$Lachesis_clusters;
my @tags2 = sort scaffold_sort keys %$RAD_tag_clusters;
print "Tags for Lachesis clusters: @tags\n";
print "Tags for RAD-tags clusters: @tags2\n";
map { die "Tags must be equal" unless $tags[$_] eq $tags2[$_] } ( 0..$#tags );


# Find the set of multiply placed scaffold in the RAD tags.
my %seen = ();
foreach my $tag (@tags) {
    my @cluster = @{$RAD_tag_clusters->{$tag}};
    map { ${$seen{$_}}{$tag}++ } @cluster;
}

my @multiply_placed = sort scaffold_sort grep { scalar keys %{$seen{$_}} > 1 } keys %seen;
print "N scaffolds multiply placed in RAD tags: ", scalar @multiply_placed, "\n";
print "@multiply_placed\n";


# Find the set of scaffolds that are placed differently in RAD tags and Lachesis tags.
my (%scaffold_to_RAD_tag_LG, %scaffold_to_Lachesis_LG);
foreach my $tag (@tags) {
    my @RAD_tags_cluster = @{$RAD_tag_clusters->{$tag}};
    my @Lachesis_cluster = @{$Lachesis_clusters->{$tag}};
    
    map { $scaffold_to_RAD_tag_LG {$_} = $tag             } @RAD_tags_cluster;
    map { $scaffold_to_Lachesis_LG{$_} = $tag unless /\D/ } @Lachesis_cluster;
}

my @inconsistently_placed;
foreach (0..$N_contigs-1) {
    if ( defined $scaffold_to_RAD_tag_LG{$_} && defined $scaffold_to_Lachesis_LG{$_} ) { # contig is placed in both sets of clusters
	if ( $scaffold_to_RAD_tag_LG{$_} != $scaffold_to_Lachesis_LG{$_} ) {
	    unless ( exists ${$seen{$_}}{$scaffold_to_Lachesis_LG{$_} } ) { # make sure contig isn't just multiply placed
		push @inconsistently_placed, $_;
		#print "INCONSISTENTLY PLACED: $_\t$scaffold_to_RAD_tag_LG{$_}\t$scaffold_to_Lachesis_LG{$_}\t",%{$seen{$_}},"\n";
	    }
	}
    }
}
print "N scaffolds inconsistently placed between RAD tags and Lachesis: ", scalar @inconsistently_placed, "\n";


# Make a "blacklist_or_multiple" lookup table: all contigs that are multiply or inconsistently placed.
my %blacklist_or_multiple = ();
map { $blacklist_or_multiple{$_}++ } ( @multiply_placed, @inconsistently_placed );
#my @inconsistently_placed = sort scaffold_sort keys %inconsistently_placed;
print "Blacklist: @inconsistently_placed\n";
print "Blacklist size: ", scalar @inconsistently_placed, "\n";



# Make a permissive merged set of scaffolds.  It will include all scaffolds placed by either Adam or me, except ones that have been multiply placed by Adam.
# Scaffolds 12, 18, 28, 20 will be reported as broken into their constituent parts.
# This file is designed for input into Lachesis via ClusterVec::ReadFile().
my $N_contigs_used = 0;
my %merged_clusters;
my (%contig_appearances);

foreach my $tag (@tags) {
    
    # Get the clusters.
    my @RAD_tags_cluster = @{$RAD_tag_clusters->{$tag}};
    my @Lachesis_cluster = @{$Lachesis_clusters->{$tag}};
    
    # Merge the clusters.  Filter out the blacklists.
    my %merged_cluster = ();
    map { $merged_cluster{$_}++ unless $blacklist_or_multiple{$_}; $contig_appearances{$_}++    } @RAD_tags_cluster;
    map { $merged_cluster{$_}++ unless $blacklist_or_multiple{$_}; $contig_appearances{$_} += 3 } @Lachesis_cluster;
    
    # Record this cluster.
    my @merged_cluster = sort scaffold_sort keys %merged_cluster;
    my $cluster_size = scalar @merged_cluster;
    #print "Cluster with tag $tag:\t$cluster_size contigs\n";
    #print "Cluster with tag $tag:\t@merged_cluster\n";
    $N_contigs_used += $cluster_size;
    $merged_clusters{$tag} = \@merged_cluster;
    
    #print "$_\t", scalar @RAD_tags_cluster, "\t", scalar @Lachesis_cluster, "\n";
}




# Print the output file!
open OUT, '>', $merged_clusters_file;

# Start with comments.
foreach (@comments) {
    s/Number of contigs in clusters\: \d+/Number of contigs in clusters\: $N_contigs_used/;
    print OUT;
}

# Now print the newly merged clusters.
foreach my $tag ( sort scaffold_sort keys %merged_clusters ) {
    my @cluster = @{$merged_clusters{$tag}};
    foreach (@cluster) {
	print OUT "\t" if ($_ ne $cluster[0]);
	#$_ = '316_and_214' if $_ eq '214'; # kludgey solution to problem of two contigs that Adam wanted to merge
	print OUT "scaffold_$_";
    }
    print OUT "\n";
}
close OUT;


# From the list of contig appearances, derive the permissive and conservative groups and the never-placed contigs.
# Being placed in Lachesis is +3 to contig appearance; each appearance in RAD tags is +1.
my (@never_seen, @permissive, @conservative, @in_Lachesis_only, @in_RAD_tags_only, @wtf );
foreach my $contig (@$scaffold_names) {
    
    if ( !exists $contig_appearances{$contig} ) { # contigs never seen
	push @never_seen, $contig;
	next;
    }
    
    push @wtf, $contig if ( exists $blacklist_or_multiple{$contig} and $blacklist_or_multiple{$contig} == 2 );
    
    next if ( $blacklist_or_multiple{$contig} ); # if a contig is on the blacklist or multiply placed, it shouldn't be in permissive or conservative groups!
    
    my $a = $contig_appearances{$contig};
    if ($a == 2 || $a == 5) { die unless $blacklist_or_multiple{$contig}; } # placed twice by RAD tags: should be on blacklist
    elsif ($a == 4) { # placed once each by Lachesis and RAD tags: ideal
	push @permissive, $contig;
	push @conservative, $contig;
    }
    elsif ( $a == 1 ) { # placed once in RAD tags, not in Lachesis
	push @in_RAD_tags_only, $contig;
	push @permissive, $contig;
    }
    elsif ( $a == 3 ) { # not placed in RAD tags, but in Lachesis
	push @in_Lachesis_only, $contig;
	push @permissive, $contig;
    }
}


print "Double-blacklist (contigs placed multiply by Adam and then discordantly by me): @wtf\n";

# Report on stuff.
report_on_scaffolds( "Entire assembly", @$scaffold_names );
report_on_scaffolds( "Scaffolds not placed in either set of linkage groups", @never_seen );
report_on_scaffolds( "Scaffolds multiply placed by Adam (i.e., scaffolds containing RAD tags in multiple linkage groups)", @multiply_placed );
report_on_scaffolds( "Scaffolds multiply placed by Adam, and inconsistently by me", @wtf );
report_on_scaffolds( "Scaffolds placed exactly once by Adam, but not by me", @in_RAD_tags_only );
report_on_scaffolds( "Scaffolds placed by me, but not by Adam", @in_Lachesis_only );
report_on_scaffolds( "'Conservative' set: scaffolds placed once by both Adam and me, and in the same group", @conservative );
report_on_scaffolds( "'Permissive' set: scaffolds placed once by Adam and/or me, with no inconsistencies", @permissive );
report_on_scaffolds( "'Blacklist': scaffolds placed once by both Adam and me, but in different groups", @inconsistently_placed );
