#!/usr/bin/env bash
# Download a GENCODE human annotation into data/. Usage: bash download_gencode.sh [version]
set -euo pipefail
VER="${1:-38}"
cd "$(dirname "$0")"
mkdir -p data && cd data
FILE="gencode.v${VER}.annotation.gtf.gz"
URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VER}/${FILE}"
if [ -f "$FILE" ]; then
  echo "$FILE already present"
else
  echo "downloading $URL"
  curl -L -o "$FILE" "$URL"
fi
echo "done: $(ls -la "$FILE")"
