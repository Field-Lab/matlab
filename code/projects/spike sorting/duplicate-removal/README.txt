---------------------------------
DUPLICATE REMOVAL README
(tamachado 2009-06-23):
---------------------------------

Start by looking at the notes about duplicate removal in the notes subfolder. They should help you to get an idea of the kinds of approaches towards duplicate finding that I was thinking about.

The first script that you should probably look at uses single linkage clustering to create duplicate groups based on an arbitrary cell comparison metric. Instead of a choosing one cell from each duplicate group, it merges all spikes in each group into a single cell. This algorithm can be run using dr_remove_duplicates_using_dendrogram. Two cell comparison metrics were implemented: one that uses spike times and another that uses EIs.

Another approach, that was more complicated and wasn't completely finished (based on a paper found in the references folder), can be found in the koosh subfolder. You can run it by executing the dr_test_framework script.