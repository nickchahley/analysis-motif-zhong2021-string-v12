#!/usr/bin/env bash
URL="https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/9606.protein.physical.links.detailed.v12.0.txt.gz"

PREFIX='data/string'
FILE=$(basename "$URL")

[[ ! -d "data/string" ]] && mkdir "$PREFIX" 
[[ -r "$PREFIX/FILE" ]] || wget "$URL" -P "$PREFIX" 
