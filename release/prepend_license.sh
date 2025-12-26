#!/bin/bash

# check for input
if [ -z "$1" ]; then
  echo "Usage: $0 path/to/dir"
  exit 1
fi

DIR="$1"

# check if license.txt exists
if [ ! -f license.txt ]; then
  echo "Error: license.txt not found in current directory"
  exit 1
fi

# prepend license.txt to each .m file in the directory tree
find "$DIR" -type f -name "*.m" | while IFS= read -r file; do
  cat license.txt "$file" > "$file.new" && mv "$file.new" "$file"
done
