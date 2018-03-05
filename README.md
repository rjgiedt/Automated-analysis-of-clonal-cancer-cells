# Automated analysis of clonal cancer cells
## Summary
Longitudinal analyses of single cell lineages over prolonged periods have been challenging particularly in processes characterized by high cell turn-over such as inflammation, proliferation, or cancer. RGB marking has emerged as an elegant approach for enabling such investigations. However, methods for automated image analysis continue to be lacking. Here, to address this, we created a number of different multicolored poly- and monoclonal cancer cell lines for in vitro and in vivo use. To classify these cells in large scale data sets, we subsequently developed and tested an automated algorithm based on hue selection. Our results showed that this method allows accurate analyses at a fraction of the computational time required by more complex color classification methods. Moreover, the methodology should be broadly applicable to both in vitro and in vivo analyses.

## Code
Included in this repository is a single matlab script, designed to take in images (placed in 3 separate channels (corresponding to microscope output rather than traditional RGB images), segment single cells/ objects from those images, and then place them into a grouping based on their color. To segment colors in this algorithm, hard cutoffs which were empirically determined were utilized to group the cells objects. 

For additional information, see https://www.ncbi.nlm.nih.gov/pubmed/24349895 to access the full manuscript.
