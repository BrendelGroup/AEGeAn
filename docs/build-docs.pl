#!/usr/bin/env perl
use strict;
use warnings;

print "AEGeAn C API\n============\n\n";

my $dir = $ARGV[0];
opendir(my $dh, $dir) or die("open dir fail");
while(my $entry = readdir($dh))
{
  if($entry =~ m/\.h$/)
  {
    process_file("$dir/$entry");
  }
}
closedir($dh);

sub process_file
{
  my $file = shift(@_);
  my $contents = do
  {
    local $/;
    open(my $fh, "<", $file) or die("open file fail: $file");
    <$fh>;
  };

  my @docblocks = $contents =~ m/(\/\*\*.+?;)/sg;
  return if(@docblocks == 0);
  
  my $descstr = shift(@docblocks);
  my ($type, $name, $desc) = parse_description($descstr);
  my $title = "$type $name";
  printf("%s\n%s\n%s\n\n", $title, "-" x length($title), $desc);

  foreach my $block(@docblocks)
  {
    parse_function_doc($block);
  }
}

sub parse_description
{
  my $desc = shift(@_);
  do {} while($desc =~ s/\s+\*\s+/ /g);
  if($desc =~ m/\@class/)
  {
    my($class_name, $description) = $desc =~ m/\/\*\* \@class (\S+) (.+)/;
    return ("Class", $class_name, $description);
  }
  elsif($desc =~ m/\@module/)
  {
    my($module_name, $description) = $desc =~ m/\/\*\* \@module (\S+) (.+)/;
    return ("Module", $module_name, $description);
  }

  die("desc fail $desc");
}

sub parse_function_doc
{
  
}
