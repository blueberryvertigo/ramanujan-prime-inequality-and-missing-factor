"""Extract text from every PDF in the directory given as argv[1] with PyMuPDF (keeps Greek letters),
writing <name>.mu.txt next to each.  Run with: python -I extract_mupdf.py pdfs"""
import sys, os, fitz
d = sys.argv[1]
for f in sorted(os.listdir(d)):
    if not f.endswith('.pdf'): continue
    t = os.path.join(d, f[:-4] + '.mu.txt')
    if os.path.exists(t) and os.path.getsize(t) > 0: continue
    try:
        doc = fitz.open(os.path.join(d, f))
        s = '\n'.join('=== page %d ===\n' % (i + 1) + p.get_text() for i, p in enumerate(doc))
        open(t, 'w', encoding='utf-8').write(s)
    except Exception as e:
        print('FAIL', f, e)
