#!/usr/bin/perl -I /booleanfs/sahoo/scripts

if (scalar(@ARGV) <= 0) {
    print "perl analyze.pl <cmd> <args>\n";
    exit(1);
}

use U;
use Hegemon;

my $cmd = shift(@ARGV);
if ($cmd eq "convert-hs") {
    &convertHS(@ARGV);
}

sub convertHS {
  my ($if, $cf, $ef, $idxf, @rest) = @_;
  my $file1 = "GSE153129_family.soft";
  my $log = 5;
  my $table1 = &U::getExprSoft($file1, undef, undef);
  my $obj = {};
  $obj->{'hash'} = {};
  $obj->{'headers'} = [];
  $obj->{'cheaders'} = ["Title"];
  $obj->{'platformids'} = [];
  $obj->{'symbols'} = {};
  $obj->{'narray'} = [];
  my $num = scalar(@{$obj->{'cheaders'}});
  $obj->{'carray'} = [0 .. ($num - 1)];

  my $cheader = 0;
  my $chh = {};
  my $ahash = {};

  foreach my $p ("GPL20301") {
  foreach my $arr (keys %{$table1->[0]->{$p}}) {
    my $l = $table1->[0]->{$p}->{$arr};
    my $list = [];
    $list->[0] = $l->[0];
    if ($cheader == 0) {
      my @data = map { if (ref($_) ne "HASH") {(split(/: /, $_))[0]} else { ""
} } grep { /: / } @$l;
      my @chdr = map { $data[$_] } 
      grep { length($data[$_]) > 1 && length($data[$_]) <= 80} 
      0 ..  $#data;
      foreach my $ch (@chdr) {
        next if $ch eq "valuefield";
        next if $ch eq "desc";
        next if defined $chh->{$ch};
        $chh->{$ch} = scalar(@{$obj->{'cheaders'}});
        push @{$obj->{'cheaders'}}, $ch;
      }
      $num = scalar(@{$obj->{'cheaders'}});
      $obj->{'carray'} = [0 .. ($num - 1)];
      for (my $i = 0; $i < $num; $i++) {
        $chh->{$obj->{'cheaders'}->[$i]} = $i;
      }
    }
  }
  }
  foreach my $p ("GPL20301") {
  foreach my $arr (keys %{$table1->[0]->{$p}}) {
    my $l = $table1->[0]->{$p}->{$arr};
    my $list = [];
    $list->[0] = $l->[0];
    $list->[1] = $p;
    foreach my $l1 (@$l) {
     foreach my $i (1 .. ($num - 1)) {
       my $str = $obj->{'cheaders'}->[$i];
       $str =~ s/[\(\)]/\./g;
       if ($l1 =~ /$str: (.*)$/) { $list->[$i] = $1; }
     }
    }
    my $title = $l->[0];
    $ahash->{$title} = $arr;
    $list->[$num] = {};
    my $time = "";
    my $status = "";
    $obj->{'hash'}->{$arr} = [$time, $status, $list];
    push @{$obj->{'headers'}}, $arr;
    print join("\t", $title, $arr), "\n";
  }
  }
  foreach my $i (0 .. ($num - 1)) {
    $obj->{'cheaders'}->[$i] =~ s/\s*\(.*\)//g;
  }

  my $file = "GSE153129_SPT6.RNAseq.txt";
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $h = <$fh>;
  $h =~ s/[\r\n]//g;
  my @headers = split("\t", $h);
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t"); 
    my $id = $list[0];
    my ($name, $desc) = split(/__/, $id);
    if (!defined $obj->{'symbols'}->{$id}) {
      $obj->{'symbols'}->{$id} = [$name, $desc];
      push @{$obj->{'platformids'}}, $id;
    }
    foreach my $i (1 .. $#headers) {
      my $aid = $headers[$i];
      if (!defined $ahash->{$aid}) {
        die "Can't find $aid\n";
      }
      my $arr = $ahash->{$aid};
      my $v = &U::getVal($list[$i], undef);
      $obj->{'hash'}->{$arr}->[2]->[$num]->{$id} = $v;
    }
  } 
  close($fh); 

  my $arr = $obj->{'headers'}->[0];
  my $l = $obj->{'hash'}->{$arr}->[2];
  my $id = $obj->{'platformids'}->[0];
  my $ex = $l->[scalar(@$l) - 1]->{$id};
  print $ex, " ", scalar(@$l), "\n";

  &U::writeObject($obj, $if, $cf, $ef, $idxf);
}

