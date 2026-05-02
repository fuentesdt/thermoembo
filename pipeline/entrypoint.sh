#!/bin/bash
# Allow the pipeline to run as any UID/GID (e.g. via --user $(id -u):$(id -g)).
# Python/VTK write cache files to $HOME; if the invoking user has no home dir
# inside the container, fall back to /tmp so those writes don't crash.
set -e
export HOME=${HOME:-/tmp}
exec python /work/pipeline/run_thermoembo.py "$@"
