# Photometry
This folder contains all the photometry data for J1816 collected.

1) AAVSO: this contains all the amateur data on J1816 (right now this consists of 20 different observers, in several bands). The data can be obtained from the AAVSO website by searching for ASASSN-V J181654.06-202117.6. The bands observed are Johnson <em>B</em> and <em>V</em>, Cousins <em>R</em> and <em>I</em>, <em>CV</em> (no filter reduced to <em>V</em>) and Sloan <em>G</em>. This contains the light curve data.

2) ASAS: this contains all the data from the ASAS survey. There is a file for the <em>V</em> band data, for the <em>I</em> band data, and the <em>V</em> band data from the newest (Nv) location. Obtained by private communication. This contains the light curve data.

3) ASAS-SN: this contains the <em>V</em> and <em>g'</em> data, which can be obtained from the survey website by searching for the coordinates of J1816. This is the survey that started this whole journey. This contains the light curve data.

4) ATLAS: this contains the <em>c</em> and <em>o</em> data from the survey. Obtained by private communication. This contains the light curve data.

5) EvryScope: this contains the data in the Sloan <em>g</em> band. Obtained by private communication. This contains the light curve data.

6) LCO-GT: this contains data in fits files that needs to be converted to a light curve. Obtained through extra time on a J0600 proposal (see https://github.com/mkenworthy/asas-sn-J060000).

7) POBS: this contains data from the Perth Observatory Research Group, provided by Craig Bowers and Rafik Bourne.

This will be compiled to one data file with one format such that only one data file is necessary and such that just one function is necessary to load the data.
