
import json
with open("/home/kelto/fungus-sv/figures/svData.json") as f: svData = json.load(f)
with open("/home/kelto/fungus-sv/figures/geneOverlap.json") as f: geneOverlap = json.load(f)
with open("/home/kelto/fungus-sv/figures/summary.json") as f: summary = json.load(f)
with open("/home/kelto/fungus-sv/figures/strainColors.json") as f: strainColors = json.load(f)
with open("/home/kelto/fungus-sv/figures/chromNames.json") as f: chromNames = json.load(f)

maxSize = 1531933
dhffc = {"total_tested":266,"confirmed":115,"weak":67,"not_del_dup":84}
inv = {"total_inv":11,"confirmed":11,"range":"32-1,091"}
cal = {"dd":0.3,"dw":0.7,"ds":2.0,"dsw":1.3,"dep":0.35,"brk":0.30,"asm":0.20,"km":0.15}
strains = ["S288C","BJ4","IMX2600","Makgeolli","SX2"]

html = open("/home/kelto/fungus-sv/figures/report_template.html").read()
html = html.replace("__SVDATA__", json.dumps(svData))
html = html.replace("__GENES__", json.dumps(geneOverlap))
html = html.replace("__SUMMARY__", json.dumps(summary))
html = html.replace("__COLORS__", json.dumps(strainColors))
html = html.replace("__CHROMS__", json.dumps(chromNames))
html = html.replace("__MAXSIZE__", str(maxSize))
html = html.replace("__STRAINS__", json.dumps(strains))
html = html.replace("__DHFFC__", json.dumps(dhffc))
html = html.replace("__INV__", json.dumps(inv))
html = html.replace("__CAL__", json.dumps(cal))

with open("/home/kelto/fungus-sv/figures/FUNGUS_SV_report.html", "w") as f:
    f.write(html)
print("Done")
