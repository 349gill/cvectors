#!/bin/bash

if [ -z "$1" ]; then
  echo "No argument supplied."
  exit 1
fi

case "$1" in
  test)
    make test
    ./test
    rm test
    ;;
  vector)
    make vector
    ./test
    rm test
    ;;
  quat)
    make quat
    ./test
    rm test
    ;;
  matrix)
    make matrix
    ./test
    rm test
    ;;
  *)
    echo "Invalid argument."
    exit 1
    ;;
esac
