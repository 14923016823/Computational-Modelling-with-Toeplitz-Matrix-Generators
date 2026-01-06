import pandas as pd
import os
import xml.etree.ElementTree as ET
import re
import matplotlib.pyplot as plt
import numpy as np
import sys

def parse_advisor_xml(path):
    if not os.path.exists(path):
        print(f"Error: Directory not found at {path}")
        return None
    
    # Helper to strip namespace prefixes (e.g., {url}name -> name)
    def clean_tag(tag):
        return re.sub(r'\{.*\}', '', tag)

    run_code = path.split('\\')[-1].strip()
    tree = ET.parse(path + rf"\{run_code}.advixe")
    root = tree.getroot()
    
    cpu_element = root.find('rdmgr').find('CPUBrandName')
    CPUmodel = cpu_element.text if cpu_element is not None else "Unknown CPU"

    tree = ET.parse(path + r"\metrics.advisum")
    root = tree.getroot()

    # 1. Parse RoofItems (Hardware Limits)
    roof_data = []
    for item in root.findall('.//roofItem'):
        # Get attributes and strip namespaces from keys
        row = {clean_tag(k): v for k, v in item.attrib.items()}
        # Get sub-element boolean types
        types = item.find('types')
        if types is not None:
            row.update({clean_tag(k): v for k, v in types.attrib.items()})
        roof_data.append(row)
    
    df_roofs = pd.DataFrame(roof_data)
    # Convert bandwidth to numeric
    if 'bandwidth' in df_roofs.columns:
        df_roofs['bandwidth'] = pd.to_numeric(df_roofs['bandwidth'])

    # 2. Parse Metrics (Application Performance)
    metrics_element = root.find('Metrics')
    if metrics_element is not None:
        metrics_data = {clean_tag(k): v for k, v in metrics_element.attrib.items()}
        df_metrics = pd.DataFrame([metrics_data])
        
        # Convert numeric columns automatically
        for col in df_metrics.columns:
            df_metrics[col] = pd.to_numeric(df_metrics[col], errors='coerce')
    else:
        df_metrics = pd.DataFrame()

    return CPUmodel, df_roofs, df_metrics

if len(sys.argv) > 1:
    path = sys.argv[1]
else:
    path = input("FULL path to the analysis directory (e000/hs000): ").strip().strip('"')

CPUmodel, df_roofs, df_metrics = parse_advisor_xml(path)

# --- 1. Setup Helper to Extract Values ---
def get_val(df, name):
    row = df.loc[df['name'] == name]
    return row['bandwidth'].values[0] if not row.empty else None

plt.figure(figsize=(16, 8))
x_ai = np.logspace(-2, 2, 500) # Range for Arithmetic Intensity

# --- 2. Plot Memory Roofs (Slanted Lines: y = x * bandwidth) ---
for roof_name in ['DRAM Bandwidth', 'L3 Bandwidth', 'L2 Bandwidth', 'L1 Bandwidth']:
    bw = get_val(df_roofs, f"{roof_name} (single-threaded)")
    if bw:
        plt.plot(x_ai, x_ai * bw, label=f'{roof_name} ({bw:.1f} GB/s)', linestyle='--')

# --- 3. Plot Compute Roofs (Horizontal Lines: y = peak) ---
peaks = {
    'DP Vector FMA Peak': 'solid',
    'DP Vector Add Peak': 'dotted',
}

for peak_name, style in peaks.items():
    peak_val = get_val(df_roofs, f"{peak_name} (single-threaded)")
    if peak_val:
        plt.axhline(y=peak_val, label=f'{peak_name} ({peak_val:.1f} GFLOPS)', 
                    linestyle=style, color='gray', alpha=0.7)

# --- 4. Plot Application Data Point ---
app_ai = df_metrics['TotalFloatAI'].values[0]
app_gflops = df_metrics['TotalGFLOPS'].values[0]

plt.scatter(app_ai, app_gflops, color='red', s=200, edgecolors='black', zorder=10)

# 5. Formatting the Log-Log Plot
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Arithmetic Intensity (FLOPs/Byte)')
plt.ylabel('Performance (GFLOPs/s)')
plt.grid(True, which="both", ls="-", alpha=0.3)
plt.legend()

plt.xlim(10**-2, 10**2)
plt.ylim(10**-1, 10**3)

plt.tight_layout()
plt.show()