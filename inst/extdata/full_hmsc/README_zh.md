# hmscHist full HMSC example data

这个文件夹用于完整跑通：

1. 读取地点 × 物种矩阵
2. 创建 hmscHist 项目对象
3. 构建和标准化历史表
4. 计算历史指数
5. 生成 HMSC 的 XData、TrData、random effects
6. 拟合 HMSC probit 模型
7. 运行短 MCMC
8. 计算预测值、模型拟合和方差分解

文件说明：

- `comm.csv`：地点 × 物种 0/1 矩阵
- `site_region.csv`：地点所属区域
- `env.csv`：现代环境变量
- `coords.csv`：经纬度
- `traits.csv`：普通物种性状
- `species_region_history.csv`：物种-区域历史表
- `speciation_events.csv`：物种形成事件表
- `history_events.csv`：历史扩散和区域丢失事件表
- `trait_history.csv`：性状历史表
- `tip_ranges.csv`：物种现代区域分布表，可选
- `phylo.tre`：toy 系统发育树

运行脚本：

```r
source(system.file("examples/full_hmsc_workflow_zh.R", package = "HmscEcoEvo"))
```

如果你使用 `devtools::load_all()` 而不是正式安装包，也可以在包根目录运行：

```r
source("inst/examples/full_hmsc_workflow_zh.R")
```

注意：示例里的 MCMC 很短，只用于确认流程可以跑通。正式分析需要更长的 MCMC。
