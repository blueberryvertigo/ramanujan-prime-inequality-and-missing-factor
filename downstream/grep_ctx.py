"""Print body lines (before the bibliography) mentioning BKLNW / Buthe / their tables, with context.
Usage: python -I grep_ctx.py pdfs ID [ID...]"""
import sys, os, re
sys.stdout.reconfigure(encoding="utf-8")
PAT = re.compile(r'Broadbent|BKLNW|Bro\+|BKL\+|B.?.?uthe|Buthe|BUTHE|Table 8\b|Table 15\b|Table 14\b|1\.93378|0\.024334|57\.184|Cor(ollary)?\.? ?11\.1|Cor(ollary)?\.? ?5\.1|Cor(ollary)?\.? ?2\.1')
for pid in sys.argv[2:]:
    lines = open(os.path.join(sys.argv[1], pid + '.mu.txt'), encoding='utf-8', errors='replace').read().split('\n')
    end = len(lines)
    for i in range(len(lines) - 1, 0, -1):
        if re.match(r'\s*(References|REFERENCES|Bibliography|BIBLIOGRAPHY)\s*$', lines[i]): end = i; break
    print('#####', pid, '(body lines', end, ')')
    last = -10
    for i in range(end):
        if PAT.search(lines[i]) and i > last + 1:
            ctx = ' '.join(x.strip() for x in lines[max(0, i - 2):i + 3])
            print(f'  l{i}: {ctx[:420]}')
            last = i
