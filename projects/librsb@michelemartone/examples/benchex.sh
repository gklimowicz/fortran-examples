#!/bin/bash

set -e
set -x

# the benchmark command assumes A.mtx is a file in Matrix Market format
rsbench -oa -Ob --bench --nmb -f pd.mtx

# it has many (librsb development-oriented) options
rsbench -oa -Ob --help

# this is mostly a development tool so don't rely on much more than the above.
