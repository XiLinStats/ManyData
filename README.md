# ManyData
This repo contains code for our article [Combining experimental and observational data through a power likelihood](https://arxiv.org/abs/2304.02339)

**Abstract**
Randomized controlled trials are commonly regarded as the gold standard for causal inference and play a pivotal role in modern evidence-based medicine. However, the sample sizes they use are often too limited to draw significant causal conclusions for subgroups that are less prevalent in the population. In contrast, observational data are becoming increasingly accessible in large volumes but can be subject to bias as a result of hidden confounding. Given these complementary features, we propose a power likelihood approach to augmenting RCTs with observational data for robust estimation of heterogeneous treatment effects. We provide a data-adaptive procedure for maximizing the Expected Log Predictive Density (ELPD) to select the influence factor that best regulates the information from the observational data. We conduct a simulation study to illustrate the efficacy of our method and its favourable features compared to existing approaches. Lastly, we apply the proposed method to data from Tennessee's Student Teacher Achievement Ratio (STAR) Study to demonstrate its usefulness and practicality in real-world data analysis.

# Installation
You can install the package ManyData from my repository:
```
devtools::install_github("XiLinStats/ManyData",build_vignettes = FALSE)
```
