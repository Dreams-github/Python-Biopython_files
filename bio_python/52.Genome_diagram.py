from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
data=SeqIO.read("NC_005816.gb","genbank")
gd_diagram=GenomeDiagram.Diagram()

gd_track_for_fet=gd_diagram.new_track(1,name="annotated features")

gd_features_set=gd_track_for_fet.new_set()

for feature in data.features:
    if feature.type!="gene":
        continue

    if len(gd_features_set)%2==0:
        color=colors.blue
    else:
        color=colors.lightblue

    gd_features_set.add_feature(feature,color=color,label=True)

gd_diagram.draw(
    format="linear",
    orientation="landscape",
    pagesize="A4",
    fragments=4,
    start=0,
    end=len(data)
)

gd_diagram.write("genome.png","PNG")

gd_diagram.draw(
    format="circular",
    circular=True,
    pagesize=(20 * cm, 20 * cm),
    start=0,
    end=len(data),
    circle_core=0.7,
)
gd_diagram.write("circular.png", "PNG")
