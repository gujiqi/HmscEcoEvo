# hmscHist 完整 HMSC 结果展示示例

这个目录里的脚本按照 HMSC 常见结果流程组织，用 toy 数据跑通：

1. `S1_define_models_hmscHist.R`：读取数据、计算 hmscHist 历史指数、生成 HMSC 输入、定义模型。
2. `S2_fit_models_short.R`：用很短的 MCMC 拟合模型。
3. `S3_evaluate_convergence.R`：输出 MCMC 收敛诊断，包括 ESS、PSRF 和 trace plots。
4. `S4_compute_model_fit.R`：输出模型拟合、交叉验证、WAIC 和预测值。
5. `S5_show_model_fit_and_variance.R`：输出方差分解、环境 vs 历史分组方差分解、物种关联和拟合图。
6. `S6_show_parameter_estimates.R`：输出 Beta、Gamma、Rho 参数估计表和热图。
7. `S7_make_predictions.R`：输出环境/历史变量梯度预测图和结果索引。

最简单运行方式：

```r
source(system.file("examples/full_hmsc_results_workflow/00_run_all_results_zh.R", package = "HmscEcoEvo"))
```

如果你在包根目录用 `devtools::load_all()`，也可以运行：

```r
source("inst/examples/full_hmsc_results_workflow/00_run_all_results_zh.R")
```

输出目录会在当前工作目录下生成：

```text
hmscHist_full_results_example/
├── data/
├── models/
└── results/
    ├── tables/
    └── plots/
```

注意：

- 示例的 MCMC 很短，只为展示所有结果文件如何生成。
- 正式分析不能使用这个 MCMC 设置解释生态结果。
- 正式分析建议提高 `samples`、`transient`、`thin` 和 `nChains`。
