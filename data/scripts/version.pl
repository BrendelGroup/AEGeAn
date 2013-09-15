#!/usr/bin/env perl
use strict;
use Getopt::Long;

sub print_usage
{
  my $OUT = shift(@_);
  my $scriptname = `basename $0`;
  chomp($scriptname);
  print $OUT "\n$scriptname: create the AgnVersion.h header file        
By default, this script reads the date and SHA1 hash of the latest commit in the
git repository. If the git repository has been removed, then values can be
manually specified.

Usage: perl $0 [options] > inc/core/AgnVersion.h
  Options:
    -d|--date: STRING       copyright date associated with this version; default
                            is the year of the latest commit in the git
                            repository; if the git repository has been removed,
                            this value must be provided explicitly
    -h|--help               print this help message and exit
    -l|--link: STRING       URL providing additional information regarding this
                            version; default is

                                https://github.com/standage/AEGeAn/commit/\$hash

                            where \$hash is the SHA1 hash of the latest commit
                            in the git repository; if the git repository has been
                            removed, this value must be provided explicitly
    -v|--version: STRING    version string associated with this version; default
                            is the SHA1 hash of the latest commit in the git
                            repository; if the git repository has been removed,
                            this value must be provided explicitly
";
}

my $version = "";
my $date    = "";
my $link    = "";

my $gitlog = `git log`;
unless($gitlog =~ m/fatal: Not a git repository/)
{
  ($version) = $gitlog =~ m/commit (\S+)/;
  ($date) = $gitlog =~ m/Date:\s+.+(\d{4}) /;
  $link = "https://github.com/standage/AEGeAn/commit/$version";
}

GetOptions
(
  "d|date=s"    => \$date,
  "h|help"      => sub { print_usage(\*STDOUT); exit(0); },
  "l|link=s"    => \$link,
  "v|version=s" => \$version,
);

if($version eq "" or $date eq "" or $link eq "")
{
  printf(STDERR "Error: no git repository detected; please provide values for ".
         "all options\n");
  print_usage(\*STDERR);
  exit(1);
}

printf("#ifndef AEGEAN_VERSION_H
#define AEGEAN_VERSION_H

#define AEGEAN_VERSION \"%s\"
#define AEGEAN_COPY_DATE \"%s\"
#define AEGEAN_LINK \"%s\"

#endif
", $version, $date, $link, $version);
