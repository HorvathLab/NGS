# Clone the scSNViz directory from the github page.

```
git init
git remote add -f origin <url>
git config core.sparseCheckout true

echo "scSNViz/" >> .git/info/sparse-checkout

git pull origin master
```

# Install using Rstudio
Open the scSNViz.Rproj file in Rstudio. Select Build>Install Package.
