import matplotlib.pyplot as plt
import numpy as np

# 1. The Dictionary (from your provided data)
advisor_results = {
    "RoofItems": [
        {"name": "DRAM Bandwidth", "value": 46.966437975, "unit": "GB/s", "multiThreaded": True, "memory": True},
        {"name": "L1 Bandwidth", "value": 4782.484872312, "unit": "GB/s", "multiThreaded": True, "memory": True},
        {"name": "L2 Bandwidth", "value": 1943.415052164, "unit": "GB/s", "multiThreaded": True, "memory": True},
        {"name": "L3 Bandwidth", "value": 475.053679116, "unit": "GB/s", "multiThreaded": True, "memory": True},
        {"name": "SP Vector FMA Peak", "value": 799.136732967, "unit": "GFLOP/s", "multiThreaded": True, "compute": True},
        {"name": "DP Vector FMA Peak", "value": 389.236965839, "unit": "GFLOP/s", "multiThreaded": True, "compute": True},
        {"name": "DP Vector Add Peak", "value": 193.124821181, "unit": "GFLOP/s", "multiThreaded": True, "compute": True},
    ],
    "Metrics": {
        "ProgramTime": 1.4134198001,
        "TotalGFLOPCount": 2.240,
    }
}

# 2. Setup Plot
plt.figure(figsize=(12, 8))
# Intensity range for the X-axis (log scale)
x_ai = np.logspace(-3, 2, 500) 

# 3. Draw the Roofs
for item in advisor_results["RoofItems"]:
    if item.get("compute") and item.get("multiThreaded"):
        # Horizontal lines for compute peaks
        plt.axhline(y=item["value"], label=item["name"], linestyle="--", alpha=0.7)
    
    if item.get("memory") and item.get("multiThreaded"):
        # Diagonal lines for memory bandwidths (y = bandwidth * x)
        y_bw = item["value"] * x_ai
        plt.plot(x_ai, y_bw, label=item["name"], linewidth=2)

# 4. Calculate Application Point (The Dot)
# N=2048, double precision (8 bytes). AI = FLOPs / Bytes
# For MVM: 2*N^2 FLOPs / 8*N^2 Bytes move = 0.25 AI
ai_point = 0.25 
# Performance = Total GFLOPs / Total Time
gflops_point = advisor_results["Metrics"]["TotalGFLOPCount"] / advisor_results["Metrics"]["ProgramTime"]

plt.scatter([ai_point], [gflops_point], color='red', s=150, zorder=10)
plt.annotate(f"  MVM Code ({gflops_point:.2f} GFLOPs/s)", (ai_point, gflops_point), 
             color='red', weight='bold', fontsize=12)

# 5. Formatting the Log-Log Plot
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Arithmetic Intensity (FLOPs/Byte)', fontsize=12)
plt.ylabel('Performance (GFLOPs/s)', fontsize=12)
plt.title('Roofline Model: MVM Performance vs. Hardware Limits', fontsize=14)
plt.grid(True, which="both", ls="-", alpha=0.3)
plt.legend(loc='lower right', frameon=True)

# Adjust axes to fit your specific CPU data
plt.xlim(0.01, 100)
plt.ylim(0.1, 10000)

plt.tight_layout()
plt.show()