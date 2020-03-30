# Principal Component Analysis Workshop

## Contact
Paul Julian | [Webpage](http://swampthingecology.org) | [Twitter](https://twitter.com/SwampThingPaul) | [Email](mailto:pauljulianphd@gmail.com)

 
## Description/contents
 - `data/`: Example data file for demonstration purposes.
 - `Presentation/`: Presentation files.
 - `./PCA_rawcode.R` : R-code as used in the presentation.
 - Other files : associated GitHub and R-project files.
 
## Preperation

To install packages necessary for this workshop run this code in your `R` console

```
pkg<-c("reshape","vegan","REdaS")
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

check.packages(pkg)

```

While not essential some of the data handling functions such as `TN_Combine()` and ploting functions such as `axis_fun()` use the `AnalystHelper` package.

```
devtools::install_github("SwampThingPaul/AnalystHelper")

```

