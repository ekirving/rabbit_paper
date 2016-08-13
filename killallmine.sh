#!/usr/bin/env bash

# kill all the processes that belong to me (except the shell)
kill $(ps aux | grep evan | grep -Pv 'sshd|bash|ps|grep' | awk '{print $2}')