import pandas as pd
import numpy as np
from pycirclize import Circos
import matplotlib.pyplot as plt
from natsort import natsorted
from pycirclize.utils import ColorCycler, load_eukaryote_example_dataset

from matplotlib import rcParams

plt.rcParams['figure.dpi'] = 300

rcParams['font.family'] = 'Arial'

translocationPath = '/IPMNPDAC_WGS/Data/'
chr_links = pd.read_csv(translocationPath + 'case7_15_translocations.csv')
sampleIDs = list(set(chr_links.sampleID))
# natural order
sampleIDs = natsorted(sampleIDs, key=lambda x: x.lstrip())
plt.rcParams["figure.figsize"] = (8,32*len(sampleIDs))

for sampleID in sampleIDs:
    # Load hg38 dataset (https://github.com/moshi4/pycirclize-data/tree/main/eukaryote/hg38)
    chr_bed_file, cytoband_file, _ = load_eukaryote_example_dataset("hg38")
    circos = Circos.initialize_from_bed(chr_bed_file, space=6)

    # Add cytoband tracks from cytoband file
    circos.add_cytoband_tracks((95, 100), cytoband_file)

    # Create chromosome color dict
    ColorCycler.set_cmap("hsv")
    chr_names = [s.name for s in circos.sectors]
    colors = ColorCycler.get_color_list(len(chr_names))
    chr_name2color = {name: color for name, color in zip(chr_names, colors)}

# Plot chromosome name & xticks
    for sector in circos.sectors:
        sector.text(sector.name, r=120, size=16, weight='bold', color=chr_name2color[sector.name])
        sector.get_track("cytoband").xticks_by_interval(40000000,label_size=12,
            label_orientation="vertical", label_formatter=lambda v: f"{v / 1000000:.0f} Mb")

    # Plot chromosome link
    chr_links_case = chr_links[chr_links.sampleID.str.contains(sampleID)][list(chr_links)[:-1]]
    caseChrLinks = chr_links_case.query('sampleID==@sampleID')[list(chr_links_case)[1:]]
    caseChrLinks['chrom1'] = ['chr'+str(a) for a in caseChrLinks.chrom1]
    caseChrLinks['chrom2'] = ['chr'+str(b) for b in caseChrLinks.chrom2]
    caseChrLinks = list(caseChrLinks.itertuples(index=False, name=None))
    for i, j in enumerate(caseChrLinks):
        region1 = j[:3]
        region2 = j[3:]
        color = chr_name2color[region1[0]]
        circos.link(region1, region2, color=color, lw=2.5)
    circos.text(sampleID, deg=315, r=150, size=16)
    fig=circos.plotfig();
