#!/usr/bin/perl
# use strict;
# use warnings;

my %clauses = ();
my $hypothesis_count = 0;
my %hypotheses = ();
my @hypParts = ();
my $thmBasis = "";
my $hypothesis = "";
my $clauseList = "";
my $proofClause = "";
my $activated = "";
my $direction = "";
my $reverse = "";
my $additionalClauseCount = 0;
my $additionalClausesString;

my %skipPredicates = ("p" => "", "e" => "");

open(my $fi,  "<",  $ARGV[0])  or die "Can't open input file: $!";
$hypothesis_count = 0;

while (<$fi>) {     # assigns each line in turn to $_
    # Get hash of clauses
    if (/^cnf\((c_[0-9]+).*?file\([^,]+,(\w+)/) { 
        $clauses{$1} = $2;
    };
    # Get array of hypotheses
    if (/\{(~\(.*)\}/) {
        $hypothesis_count++;
        $hypothesis = sprintf "hypothesis_%03d", $hypothesis_count ;
        @hypParts = ();
        $additionalClausesString = "";
        $additionalClauseCount = 0;
        foreach my $val (split(';', $1)){
            $additionalClauseCount++;
            $val =~ s/~\((.*)\)/~$1/;
            my $component = $val =~ s/^~//r;
            push(@hypParts, $val);
            # $additionalClausesString .= sprintf "cnf(%s_component_%03d,axiom,%s).\n", $hypothesis, $additionalClauseCount, $component;
            $additionalClausesString .= sprintf "%d cnf(%s_component_%03d,axiom,%s).\n", $hypothesis_count, $hypothesis, $additionalClauseCount, $component;
        };
        # $hypothesis = "(" . join("|", @hypParts) . ")";
        # print($hypothesis, "\n  ");
        # $additionalClausesString .= "\n";
        $hypotheses{$hypothesis} = $additionalClausesString;
    };
}
for(keys %hypotheses){
	# print($_);
    print($hypotheses{$_});
}

# seek $fi, 0, 0;
# while (<$fi>) {     # assigns each line in turn to $_
#     if (/^cnf\((c_[0-9]+),plain,(.*?),inference\(prop_impl_just,\[status\(thm\)\],\[([c_0-9,]+)\]/) {
#         $proofClause = $1;
#         $thmBasis = $2;
#         $clauseList = $3;
#         if (exists($hypotheses{$thmBasis})) {
#             print( '"', join('","', $thmBasis , $proofClause, $clauseList), '"', "\n");
#             print("Exists");
#             foreach my $clause (split(",", $clauseList)) {
#                 if (exists($clauses{$clause})) {
#                     $activated = $clauses{$clause};
#                     $direction = "0";

#                     if (exists $skipPredicates{substr($activated, 0, 1)}) {
#                         next;
#                     } elsif (substr($activated, 0, 1) eq "r") {
#                         if (substr($activated, length($activated) - 3, 3) ne "_in") {
#                             next;
#                         }
#                         if ( substr($activated, length($activated) - 10, 7) eq "reverse" ) {
#                             $direction = "-1";
#                         } else {
#                             $direction = "1";
#                         }

#                     }

#                     print( '"', join('","', $thmBasis , $proofClause, $activated, $direction), '"', "\n");
#                 }
#             };
#             delete($hypotheses{$thmBasis});
#         }
#     }
# }

# close $fi or die "$fi: $!";
# # close $fo or die "$fo: $!";