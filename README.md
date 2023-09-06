# machuruku
Machuruku is an R package for reconstructing the ancestral distributions of lineages with phylogenetic niche modeling. The package takes present-day occurrence data, a time-calibrated phylogeny, and past and present climate data to infer and visualize the ancient niches of species. 
## September 6 2023: Machuruku 2.0.1 released
Small patch updating `machu.tree.unc` to account for more variation in format of a multi-tree nexus file.
### August 7 2023: Machuruku 2.0 released
The new version of Machuruku is now available via Github. Machuruku 2.0 has been rewritten from the ground up, with quicker code, cleaner graphics, and the ability to reconstruct and project multiple timeslices and paleoclimates all in one command. Visit the revamped [tutorial](https://github.com/wxguillo/machuruku/tree/main/tutorial#machuruku-the-tutorial-20) for an in-depth guide at using Machuruku 2.0.
### Install
To install Machuruku, simply use the `install_github()` function from `devtools`:
```
install.packages("devtools")
devtools::install_github("wxguillo/machuruku")
```
You may encounter an error trying to install the "Treeio" package, which is not on CRAN. Install it with the following and try the above installation again:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("treeio")
```
Load up Machuruku by running:
```
library(machuruku)
```
You can test if the function worked by typing `machu` into the console (in RStudio) and seeing if all of the functions appear in autofill.

### *Software citation and associated manuscript:*
>  Guillory WG, Brown JL. 2021. A new method for integrating ecological niche modeling with phylogenetics to estimate ancestral distributions. Systematic Biology, 70(5):1033-1045. [Link.](https://academic.oup.com/sysbio/advance-article-abstract/doi/10.1093/sysbio/syab016/6171196) (Email me at wxg1@rutgers.edu for PDF)

![frog logo](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/machurukuLogoShamelessFrog.jpg)
