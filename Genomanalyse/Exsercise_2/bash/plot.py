#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

daten = ('Human Intron', 'Maus Introns', 'Human Exons', 'Maus Exon',
         'Human IntergenicRegion', 'Maus IntergenicRegion')
anzahl = [37.7, 34.22, 6.98, 5.48, 55.23, 60.28]

y_pos = np.arange(len(daten))

plt.bar(y_pos, anzahl, align='center')
plt.xticks(y_pos, daten)
plt.ylabel('Prozent')
plt.ylim(0, 100)
plt.title('Vergleich Human & Maus')
# Striche auf x-Achse ausschalten
plt.tick_params(
    axis='x',
    which='both',  # major und minor ticks
    bottom=False  # ticks auf der x-Achse (unten)
)

plt.show()
