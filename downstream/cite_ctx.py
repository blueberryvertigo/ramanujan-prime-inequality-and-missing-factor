"""For each extracted paper text pdfs/<ID>.mu.txt, find the bibliography labels of BKLNW and of
Buthe 2016 (partial RH) and print every body line citing them, with one line of context.
Usage: python -I cite_ctx.py pdfs ID [ID...]"""
import sys, os, re
sys.stdout.reconfigure(encoding="utf-8")
d = sys.argv[1]
LAB = re.compile(r'^\s*\[([^\]]{1,20})\]\s')
for pid in sys.argv[2:]:
    t = open(os.path.join(d, pid + '.mu.txt'), encoding='utf-8', errors='replace').read()
    lines = t.split('\n')
    labels = {}
    for i, L in enumerate(lines):
        m = LAB.match(L)
        if not m: continue
        blk = L
        j = i + 1
        while j < len(lines) and j < i + 4 and not LAB.match(lines[j]):
            blk += ' ' + lines[j]; j += 1
        if re.search(r'Broadbent|BROADBENT', blk): labels[m.group(1)] = 'BKLNW'
        elif re.search(r'(B.{0,3}uthe|BUTHE|B.{0,3}UTHE|Buthe)', blk) and re.search(r'partial|PARTIAL|Partial', blk):
            labels[m.group(1)] = 'Buthe16'
    print('#####', pid, labels)
    for lab, name in labels.items():
        pat = re.compile(r'\[(?:[^\]]*[ ,;])?' + re.escape(lab) + r'(?:[ ,;][^\]]*)?\]')
        for i, L in enumerate(lines):
            if pat.search(L) and not LAB.match(L):
                ctx = ' '.join(x.strip() for x in lines[max(0, i - 1):i + 2])
                print(f'  {name} l{i}: {ctx[:350]}')
