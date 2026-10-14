#!/bin/sh
# usage: fetch_arxiv.sh OUTDIR id1 id2 ...   (downloads arXiv PDFs, extracts text with pdftotext)
out="$1"; shift
mkdir -p "$out"
for id in "$@"; do
  f="$out/$(echo $id | tr '/' '_')"
  if [ ! -s "$f.pdf" ]; then
    curl -sL -A "Mozilla/5.0 (research script)" "https://arxiv.org/pdf/$id" -o "$f.pdf"
    sleep 3
  fi
  pdftotext -layout "$f.pdf" "$f.txt" 2>/dev/null || echo "FAIL $id"
done
