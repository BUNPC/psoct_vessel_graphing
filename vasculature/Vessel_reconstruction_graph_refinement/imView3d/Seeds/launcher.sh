#!/bin/sh
until test -f quitjob; do
   if [ -f dothings ]; then
        echo "roger that!"
        rm -f dothings
        # command to run
	qqsub.pl launchScript2.bat 15
        touch done
   fi
   sleep 5
   echo "check"
done

rm -f quitjob dothings
