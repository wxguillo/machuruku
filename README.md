# machuruku
## Sept 16 2021: machuruku has been updated to v 1.8.3, please reinstall.
Welcome to our page for Machuruku, our R package for phylogenetic niche modeling. Machuruku uses a modified Bioclim niche-modeling method and ancestral character estimation to reconstruct ancestral niches. Input data consists of present-day climate rasters, a time-calibrated phylogeny, and taxon occurrence data for all tips. Ancestral niche models can be projected into additionally provided paleoclimate data to visualize the geographic origins of a lineage. 
To install Machuruku, simply run the following in R:
```
install.packages("devtools")
devtools::install_github("wxguillo/machuruku")
```
If your install fails because of ggtree, install it separately with the following, then reattempt the above:
```
install.packages("BiocManager")
BiocManager::install("ggtree")
```

For both quick and detailed tutorials on how to use Machuruku, please visit the [tutorial page](https://github.com/wxguillo/machuruku/tree/main/tutorial) in this repository.

*Software citation and associated manuscript:* Guillory WG, Brown JL. 2021. A new method for integrating ecological niche modeling with phylogenetics to estimate ancestral distributions. Systematic Biology, 70(5):1033-1045. [Link.](https://academic.oup.com/sysbio/advance-article-abstract/doi/10.1093/sysbio/syab016/6171196) (Email me at wxg1@rutgers.edu for PDF)
![frog logo](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/machurukuLogoShamelessFrog.jpg)
