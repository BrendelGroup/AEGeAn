#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# Daniel S. Standage
# 13 Sept 2013
# daniel.standage@gmail.com
#
# Usage: perl build-docs.pl inc/core inc/ParsEval > docs/api.rst

print "AEGeAn C API\n============\n\nThe AEGeAn Toolkit relies heavily on data
types implemented by the GenomeTools library. For data types beginning with
``Gt``, see the GenomeTools API documentation at
http://genometools.org/libgenometools.html.\n\n";

# Main procedure
my $dir;
foreach my $dirname(@ARGV)
{
  $dir = $dirname;
  opendir(my $dh, $dir) or die("open dir fail");
  while(my $entry = readdir($dh))
  {
    make_docs_from_file("$dir/$entry") if($entry =~ m/\.h$/);
  }
  closedir($dh);
}

# Generate class or module documentation from a header file
sub make_docs_from_file
{
  my $file = shift(@_);
  my $contents = do
  {
    local $/;
    open(my $fh, "<", $file) or die("open file fail: $file");
    <$fh>;
  };

  my @docblocks = $contents =~ m/(\/\*\*.+?\*\/.+?;)/sg;
  return if(@docblocks == 0);

  process_doc_blocks(@docblocks);
}

# Process a set of doc blocks from a class/module header file
sub process_doc_blocks
{
  my $summary = shift(@_);
  do {} while($summary =~ s/\s+\*\s+/ /g);

  if($summary =~ m/\@class/)
  {
    my($class, $description) = $summary =~ m/\/\*\* \@class (\S+) (.+)/;
    my $title = "Class $class";
    my $url = "https://github.com/standage/AEGeAn/blob/master/$dir/$class.h";
    printf("%s\n%s\n\n", $title, "-" x length($title));
    printf(".. c:type:: %s\n\n  %s See the `%s class header <%s>`_.\n\n",
           $class, $description, $class, $url);

    foreach my $block(@_)
    {
      process_doc_block($block);
    }
  }

  elsif($summary =~ m/\@module/)
  {
    my($module, $description) = $summary =~ m/\/\*\* \@module (\S+) (.+)/;
    my $title = "Module $module";
    my $url = "https://github.com/standage/AEGeAn/blob/master/$dir/$module.h";
    printf("%s\n%s\n\n", $title, "-" x length($title));
    printf("%s See the `%s module header <%s>`_.\n\n", $description, $module,
           $url);

    foreach my $block(@_)
    {
      process_doc_block($block);
    }
  }

  else
  {
    die("class/module fail");
  }
}

# Process the doc block corresponding to an individual function or type def
sub process_doc_block
{
  my $block = shift(@_);

  $block =~ s/\/\*\*\W+//s;
  my($type) = $block =~ m/(function|functype|type)/;
  die("block type fail") unless($type);
  $block =~ s/(function|functype|type) *//;
  my($synopsis) = $block =~ m/^(.+?)(\@param|\@member|\*\/)/s;
  $synopsis =~ s/\*//g;
  $synopsis =~ s/\s+/ /sg;
  $synopsis =~ s/\s+$//;

  if($type eq "function")
  {
    my($prototype) = $block =~ m/\*\/\s*(\w[^;]+);/;
    $prototype =~ s/\s+/ /sg;
    print(".. c:function:: $prototype\n\n  $synopsis\n\n");
  }
  elsif($type eq "type")
  {
    my($typetype, $typename) = $block =~ m/(struct|enum) (Agn\S+)/;
    print(".. c:type:: $typename\n\n  $synopsis\n\n");
    my @members = $block =~ m/\@member (.+)/g;
    foreach my $member(@members)
    {
      my($mtype, $mname, $mdesc) = $member =~ m/\[(.+)\] (\S+) (.+)/;
      print("  * **$mtype $mname**: $mdesc\n");
    }
    print("\n\n");
  }
  elsif($type eq "functype")
  {
    my($signature) = $block =~ m/\*\/\s*(\w[^;]+);/;
    $signature =~ s/\s+/ /sg;
    print(".. c:type:: $signature\n\n  $synopsis\n\n");
  }
  else
  {
    die("block type fail");
  }
}
