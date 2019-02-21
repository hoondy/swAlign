#!/usr/bin/perl

### Author:     Donghoon Lee (donghoon (dot) lee (at) yale (dot) edu)
### Created:    02/11/2015
### Usage:      perl runSW.pl <input file> <score file>
### Example:    perl runSW.pl input.txt blosum62.txt
### Note:       Smith-Waterman Algorithm

use strict;
use warnings;
use 5.010; # minimum perl version required
use Data::Dumper qw(Dumper); # used for debugging

# usage statement
die "usage: $0 <Input Sequences> <Scoring Matrix>\n" unless @ARGV == 2;

# open input files
open my $fh1, "<$ARGV[0]" or die "Could not open $ARGV[0]: $!"; # input sequences file
open my $fh2, "<$ARGV[1]" or die "Could not open $ARGV[1]: $!"; # scoring matrix file

# read sequences from input file
my @input_seq;
while(<$fh1>){
    $_ =~ s/\r//; # removes carriage-return if there's any (i.e., sample-input.txt had them)
    chomp; 
    push @input_seq, $_;
} 
close $fh1;
my ($seq1, $seq2)=@input_seq;

# print input sequences
print "*********************\n";
print "*  Input Sequences  *\n";
print "*********************\n";
print "\n";
print "Seq 1: $seq1\n";
print "Seq 2: $seq2\n";
print "\n";

# read scoring matrix (i.e., blosum62.txt file)
my @score_mat;
while(<$fh2>){ ## read single line from the file
    chomp; ## remove newline
    push @score_mat, $_;
} 
close $fh2;

# create match/mismatch score look-up array
my %charH;
my %charV;
my @score_table;
for my $row (0 .. $#score_mat){ ## for each row of score matrix..
    my @fields=split("\\s+",$score_mat[$row]);
    for my $col (0 .. $#fields){ ## for each column of score matrix..
        if($row==0){
            $charH{$fields[$col]}=$col; ## first row, make look up table to get column
        }
        else{
            if($col==0){
                $charV{$fields[$col]}=$row; ## first column, look up table to get row
            }
            else{
                $score_table[$row][$col]=$fields[$col];
            }
        }
    }
}
# print "$charH{A}\n";
# print $score_table[$charH{A}][$charV{A}]."\n";
# print Dumper \@score_table;

# display score matrix
#print join "\n", @score_mat;

# initialize H matrix[row][column] by setting the first row and column with 0s
my @matrix;
for(my $j=0; $j<=length($seq1); $j++){
    $matrix[0][$j]{score}=0; # tracks similarity score
    $matrix[0][$j]{trace}="NA"; # contains trace direction
    $matrix[0][$j]{del}=0; # tracks previous deletion (1 = prev gap)
    $matrix[0][$j]{ins}=0; # tracks previous insertion (1 = prev gap)
    $matrix[0][$j]{match}=0; # tracks match status (1 = match)
}
for(my $i=0; $i<=length($seq2); $i++){
    $matrix[$i][0]{score}=0;
    $matrix[$i][0]{trace}="NA";
    $matrix[$i][0]{del}=0;
    $matrix[$i][0]{ins}=0;
    $matrix[$i][0]{match}=0;
}

# evaluate similarity matrix
my ($max_i, $max_j, $max_score)=(0, 0, 0);
for(my $i=1; $i<=length($seq2); $i++){ # for each row i
    for(my $j=1; $j<=length($seq1); $j++){ # for each column j
        my ($diag_score, $ins_score, $del_score);
        my $letterH=substr($seq1, $j-1, 1); # get letter from seq1 (Horizontal)
        my $letterV=substr($seq2, $i-1, 1); # get letter from seq2 (Vertical)
        
        # calculate match/mismatch score
        #print "$letterH x $letterV : ";
        #print $score_table[$charH{$letterH}][$charV{$letterV}]."\n"; # look up score
        $diag_score=$matrix[$i-1][$j-1]{score} + $score_table[$charH{$letterH}][$charV{$letterV}];
        
        # calculate gap penalties: opening gap -2, extension gap -1 
        if($matrix[$i-1][$j]{del} == 0){
            $del_score=$matrix[$i-1][$j]{score} - 2; # deletion gap, opening
        }
        else{
            $del_score=$matrix[$i-1][$j]{score} - 1; # deletion gap, extension
        }
        if ($matrix[$i][$j-1]{ins} == 0){
            $ins_score=$matrix[$i][$j-1]{score} - 2; # insertion gap, opening
        }
        else{
            $ins_score=$matrix[$i][$j-1]{score} - 1; # insertion gap, extension
        }
        
        # calculate max score
        if($diag_score<=0 and $del_score<=0 and $ins_score<=0){ # set to minimum 0 if all scores are negative
            $matrix[$i][$j]{score}=0;
            $matrix[$i][$j]{trace}="NA";
            $matrix[$i][$j]{del}=0;
            $matrix[$i][$j]{ins}=0;
            $matrix[$i][$j]{match}=0;
        }
        elsif($diag_score>=$del_score and $diag_score>=$ins_score){
            $matrix[$i][$j]{score}=$diag_score;
            $matrix[$i][$j]{trace}="diagonal";
            
            if($diag_score==$del_score){ # extend del gap if scores are the same
                $matrix[$i][$j]{del}=1;
            }
            else{
                $matrix[$i][$j]{del}=0;
            }

            if($diag_score==$ins_score){ # extend ins gap if scores are the same
                $matrix[$i][$j]{ins}=1;
            }
            else{
                $matrix[$i][$j]{ins}=0;
            }

            if($letterH eq $letterV){
                $matrix[$i][$j]{match}=1; # match
            }
            else{
                $matrix[$i][$j]{match}=0; # mismatch
            }
        }
        elsif($del_score>=$ins_score){
            $matrix[$i][$j]{score}=$del_score;
            $matrix[$i][$j]{trace}="up";
            $matrix[$i][$j]{del}=1;
            $matrix[$i][$j]{ins}=0;
            $matrix[$i][$j]{match}=0;
        }
        else{
            $matrix[$i][$j]{score}=$ins_score;
            $matrix[$i][$j]{trace}="left";
            $matrix[$i][$j]{del}=0;
            $matrix[$i][$j]{ins}=1;
            $matrix[$i][$j]{match}=0;
        }

        # record max score for trace
        if ($matrix[$i][$j]{score}>$max_score){
            $max_i=$i;
            $max_j=$j;
            $max_score=$matrix[$i][$j]{score};
        }
    }
}

# print similarity matrix
#print Dumper \@matrix;
print "******************\n";
print "*  Score Matrix  *\n";
print "******************\n";
print "\n";
for(my $j=0; $j<=length($seq1); $j++){ # horizontal
    if($j == 0){
        print "\t";
    }
    else{
        print "\t".substr($seq1, $j-1, 1);
    }
}
print "\n";
for (my $i=0; $i<=length($seq2); $i++){ # vertical
    if($i == 0){
        print"\t";
    }
    else{
        print substr($seq2, $i-1, 1)."\t";
    }
    
    for(my $j=0; $j <= length($seq1); $j++){
        #print "$i\t$j\n";
        print "$matrix[$i][$j]{score}\t";
    }
    print "\n";
}
print "\n";

# traceback
my $align1="";
my $align2="";
my $between="";
my $j=$max_j;
my $i=$max_i;
while($matrix[$i][$j]{trace} ne "NA"){
    if($matrix[$i][$j]{trace} eq "diagonal"){
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= substr($seq2, $i-1, 1);
        if($matrix[$i][$j]{match} == 1){
            $between .= "|";    
        }
        else{
            $between .= " ";
        }
        $i--;
        $j--;
    }
    elsif($matrix[$i][$j]{trace} eq "up"){ # deletion
        $align1 .= "-";
        $align2 .= substr($seq2, $i-1, 1);
        $between .= " ";
        $i--;
    }
    elsif($matrix[$i][$j]{trace} eq "left"){ # insertion
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= "-";
        $between .= " ";
        $j--;
    }   
}

# process local alignments
$align1=reverse $align1; # after tracing back, reverse the sequence
$align2=reverse $align2; # after tracing back, reverse the sequence
$between=reverse $between; # after tracing back, reverse the sequence
my $left_align1=substr($seq1, 0, $j);
my $left_align2=substr($seq2, 0, $i);
my $right_align1=substr($seq1, $max_j, length($seq1));
my $right_align2=substr($seq2, $max_i, length($seq2));
if(length($left_align1) > length($left_align2)){ # pad with spaces if needed
    $align1=$left_align1."(".$align1.")".$right_align1;
    $align2=(" " x (length($left_align1) - length($left_align2))).$left_align2."(".$align2.")".$right_align2;
    $between=(" " x (length($left_align1)))." ".$between;
}
elsif(length($left_align1) < length($left_align2)){
    $align1=(" " x (length($left_align2) - length($left_align1))).$left_align1."(".$align1.")".$right_align1;
    $align2=$left_align2."(".$align2.")".$right_align2;
    $between=(" " x (length($left_align2)))." ".$between;
}
else{
    $align1=$left_align1."(".$align1.")".$right_align1;
    $align2=$left_align2."(".$align2.")".$right_align2;
    $between=(" " x (length($left_align2)))." ".$between;   
}

# display best local alignment
print "**************************\n";
print "*  Best Local Alignment  *\n";
print "**************************\n";
print "\n";
print "Alignment Score:\t$max_score\n";
while(length($align1) > 0 or length($align2) > 0 or length($between) > 0){ # display 100 characters per line
    if(length($align1)>100){
        print substr($align1, 0, 100)."\n";
        $align1=substr($align1, 100);
    }
    elsif(length($align1) > 0){
        print $align1."\n";
        $align1="";
    }
    else{
        print "\n";
    }
    if(length($between)>100){
        print substr($between, 0, 100)."\n";
        $between=substr($between, 100);
    }
    elsif(length($between) > 0){
        print $between."\n";
        $between="";
    }
    else{
        print "\n";
    }
    if(length($align2)>100){
        print substr($align2, 0, 100)."\n";
        $align2=substr($align2, 100);
    }
    elsif(length($align2) > 0){
        print $align2."\n";
        $align2="";
    }
    else{
        print "\n";
    }
}
