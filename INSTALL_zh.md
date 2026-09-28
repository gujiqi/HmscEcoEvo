# HmscEcoEvo 本地安装说明

从源码目录安装：

```r
remotes::install_local("HmscEcoEvo", upgrade = "never", dependencies = TRUE)
library(HmscEcoEvo)
```

开发模式加载：

```r
devtools::load_all("HmscEcoEvo")
```

如果没有安装 `devtools`，也可以用 base R：

```r
R CMD INSTALL HmscEcoEvo
```

`hmscHist_data()` 等旧函数名仍然可用，但推荐新项目对象入口：

```r
proj <- HmscEcoEvo_data(comm, site_region)
```
