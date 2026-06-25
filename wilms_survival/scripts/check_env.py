"""
Environment check — confirms all required packages import correctly.
Run from the project root:  python wilms_survival/scripts/check_env.py
"""

import sys

packages = [
    ("requests",    "requests",    "__version__"),
    ("pandas",      "pandas",      "__version__"),
    ("numpy",       "numpy",       "__version__"),
    ("matplotlib",  "matplotlib",  "__version__"),
    ("lifelines",   "lifelines",   "__version__"),
    ("scipy",       "scipy",       "__version__"),
    ("seaborn",     "seaborn",     "__version__"),
]

print(f"Python  {sys.version.split()[0]}")
print("-" * 36)

all_ok = True
for display, module, attr in packages:
    try:
        mod = __import__(module)
        version = getattr(mod, attr, "?")
        print(f"  {display:<12} {version}")
    except ImportError as e:
        print(f"  {display:<12} MISSING — {e}")
        all_ok = False

print("-" * 36)
print("All packages OK" if all_ok else "SOME PACKAGES MISSING — re-run pip install")
