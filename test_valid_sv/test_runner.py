#!/usr/bin/env python3
"""Minimal test to verify VALID-SV script structure"""

import sys
import os

# Test 1: Check if the main script exists
script_path = "valid_sv/run_validation.py"
if os.path.exists(script_path):
    print(f"✓ Found {script_path}")
else:
    print(f"✗ {script_path} not found")
    sys.exit(1)

# Test 2: Syntax check
import subprocess
result = subprocess.run(["python3", "-m", "py_compile", script_path], 
                       capture_output=True, text=True)
if result.returncode == 0:
    print("✓ Syntax check passed")
else:
    print("✗ Syntax errors found:")
    print(result.stderr)
    sys.exit(1)

# Test 3: Check import statements (without executing)
print("\nChecking import structure...")
with open(script_path) as f:
    content = f.read()
    imports = [line.strip() for line in content.split('\n') 
               if line.startswith('from valid_sv') or line.startswith('import valid_sv')]
    for imp in imports:
        print(f"  {imp}")

print("\n✓ Script structure looks valid")
print("\nNote: Full execution will fail without:")
print("  - valid_sv module installed")
print("  - Input files (VCF, BAM, reference)")
print("  - Dependencies (pysam, numpy, etc.)")
