# Clone the scSNViz directory from the github page.

```
git init <repo>
cd <repo>
git remote add -f origin <url>

git config core.sparseCheckout true

# This tells git which directories you want to checkout. Then you can pull just those directories
echo "scSNViz/" >> .git/info/sparse-checkout

git pull origin master
```

# Install using Rstudio
Open the scSNViz.Rproj file in Rstudio. Select Build>Install Package.
