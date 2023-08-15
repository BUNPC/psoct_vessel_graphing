#!/usr/pubsw/bin/perl

use strict;

my $arg=$ARGV[0];
my $maxnum=$ARGV[1];
my $user=$ENV{'USER'};
my ($jobnum,$i);

if($arg eq "-h" || $arg eq "--help")
{
	print <<EOF;
	qqsub.pl -- a perl script for submitting limited number of parallel jobs

	format:
		qqsub.pl commandfile maxjobnum
	where 
	"commandfile" is a text file, each line is a command to run
	"maxjobnum"   is an integer specifing the maximum parallel jobs

	example:
		qqsub.pl myjoblist 10
EOF
	exit;
}
if($maxnum==0) {$maxnum=15;}

$jobnum=0;
$i=0;

open(JOBFILE,"<$arg") || die("can not open file");

while(<JOBFILE>)
{
	print $_." $jobnum $user\n";
	chop;
	if(length($_)==0 || not -f $_ || not ($_ =~ /^#/)) {next;}
	print "submitting ".$_."\n";
	print "pbsubmit -c \"".$_." > log/cmd".$i."_$arg.log\"\n";
	system("pbsubmit -l nodes=1:opteron -c \"".$_." > log/cmd".$i."_$arg.log\"");

        sleep 20;

        $jobnum = `qstat | grep $user | wc -l`;
        $jobnum=~s/[\r\n ]//g;

	while($jobnum>=$maxnum)
	{
	        sleep 10;
		$jobnum = `qstat | grep $user | wc -l`;
		$jobnum=~s/[\r\n ]//g;
		print "$jobnum job(s) running\n";
	}
	$i++;
}
while($jobnum>0)
{
        sleep 10;
        $jobnum = `qstat | grep $user | wc -l`;
        $jobnum=~s/[\r\n ]//g;
        print "waiting for $jobnum job(s) to complete...\n";
}

close(JOBFILE);

