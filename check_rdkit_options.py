"""Check all available RDKit drawing options"""
from rdkit.Chem.Draw import rdMolDraw2D

def list_rdkit_options():
    """List all available drawing options in RDKit"""
    drawer = rdMolDraw2D.MolDraw2DSVG(300, 250)
    opts = drawer.drawOptions()
    
    print("🎨 RDKit Drawing Options\n")
    print("="*60)
    
    # Get all attributes
    all_attrs = dir(opts)
    
    # Filter out private/magic methods
    public_attrs = [attr for attr in all_attrs if not attr.startswith('_')]
    
    # Categorize
    aromatic_related = []
    carbon_related = []
    other = []
    
    for attr in public_attrs:
        attr_lower = attr.lower()
        if 'aromatic' in attr_lower or 'circle' in attr_lower or 'ring' in attr_lower:
            aromatic_related.append(attr)
        elif 'carbon' in attr_lower or 'methyl' in attr_lower or 'atom' in attr_lower:
            carbon_related.append(attr)
        else:
            other.append(attr)
    
    print("\n🔵 AROMATIC/CIRCLE RELATED:")
    for attr in aromatic_related:
        try:
            value = getattr(opts, attr)
            print(f"   • {attr}: {value} ({type(value).__name__})")
        except:
            print(f"   • {attr}: (method)")
    
    print("\n🟢 CARBON/ATOM RELATED:")
    for attr in carbon_related:
        try:
            value = getattr(opts, attr)
            print(f"   • {attr}: {value} ({type(value).__name__})")
        except:
            print(f"   • {attr}: (method)")
    
    print("\n⚪ OTHER OPTIONS:")
    for attr in other[:20]:  # Show first 20
        try:
            value = getattr(opts, attr)
            print(f"   • {attr}: {value} ({type(value).__name__})")
        except:
            print(f"   • {attr}: (method)")
    
    print(f"\n   ... and {len(other) - 20} more options")
    
    print("\n" + "="*60)
    print("\n🔍 KEY FINDINGS:")
    
    # Check specific features
    has_aromatic_circles = any('circle' in attr.lower() for attr in aromatic_related)
    has_show_carbons = any('carbon' in attr.lower() or 'methyl' in attr.lower() for attr in carbon_related)
    
    print(f"   Aromatic Circles: {'✅ YES' if has_aromatic_circles else '❌ NO'}")
    print(f"   Show Carbons: {'✅ YES' if has_show_carbons else '❌ NO'}")

if __name__ == "__main__":
    list_rdkit_options()
