#!/bin/bash
set -eo pipefail

python -u /opt/files/driver.py | tee /data/logs/chromoseq.log
