# Install scSNViz from GitHub

### Option 1:

```
library(devtools)
install_github("HorvathLab/NGS", ref = "feature/scSNViz_R_v1.0.0", subdir = "scSNViz")
```
If the above fails due to rate limits, try generating a GitHub Personal Access Token (PAT), add it into your environment and then run again. 

Another way to do this is to configure R to use the Windows Internet API for download: 

### Option 2:
```
options(download.file.method = "wininet")   # can try other methods such as 'libcurl', 'wget', etc.
install_github("HorvathLab/NGS", ref = "feature/scSNViz_R_v1.0.0", subdir = "scSNViz")
```

### Option 3:
```
git init
git remote add -f origin <url>
git config core.sparseCheckout true

echo "scSNViz/" >> .git/info/sparse-checkout

git pull origin feature/scSNViz_R_v1.0.0
```

# Install package using Rstudio
Open the scSNViz.Rproj file in Rstudio. Select Build>Install Package.
