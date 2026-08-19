import pandas as pd, glob, matplotlib.pyplot as plt, numpy as np, os, jinja2
os.makedirs("results", exist_ok=True)

def load():
    rows=[]
    for f in glob.glob("results/*.csv"):
        df=pd.read_csv(f)
        rows.append(df)
    return pd.concat(rows, ignore_index=True)

df=load().sort_values("time_ns")
print(df)
plt.figure(figsize=(8,4)); 
plt.bar(df['name'], df['time_ns']); plt.ylabel('ns/call'); plt.xticks(rotation=45, ha='right'); plt.tight_layout()
plt.savefig("results/time.png", dpi=160)

plt.figure(figsize=(8,4)); 
plt.bar(df['name'], df['ulp95']); plt.ylabel('ULP@95%'); plt.xticks(rotation=45, ha='right'); plt.tight_layout()
plt.savefig("results/ulp95.png", dpi=160)

tpl = """<!doctype html><html><head><meta charset=\"utf-8\"><title>Benchmark log/raíz</title>
<style>body{font-family:sans-serif} table{border-collapse:collapse} td,th{border:1px solid #ccc;padding:6px}</style>
</head><body>
<h1>Benchmark log/raíz</h1>
<p>ns/call y ULP@95% (mismo hardware, software distinto).</p>
<h2>Resultados</h2>
{{ table }}
<h2>Tiempo (ns/call)</h2><img src=\"results/time.png\">
<h2>Exactitud (ULP@95%)</h2><img src=\"results/ulp95.png\">
</body></html>"""
html = jinja2.Template(tpl).render(table=df.to_html(index=False))
open("report.html","w").write(html)
print("OK -> report.html")
