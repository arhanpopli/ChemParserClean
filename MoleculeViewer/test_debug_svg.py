import sys
sys.path.insert(0, '.')
from app.chemistry import smiles_to_svg
import re

error, svg = smiles_to_svg('CCC', 400, 300, {'show_methyls': True})

with open('test_ccc_debug.svg', 'w', encoding='utf-8') as f:
    f.write(svg)
    
print('Saved test_ccc_debug.svg')
texts = re.findall(r'<text[^>]*>([^<]+)</text>', svg)
print(f'Text elements: {texts}')
print(f'Has CH₃: {"CH₃" in svg or "CH3" in svg}')
print(f'SVG snippet with CH:')
if 'CH' in svg:
    idx = svg.index('CH')
    print(svg[max(0,idx-100):idx+100])
