# Methods in R

R pipelines for the scRNA pipelines benchmark, one omnibenchmark (api 0.8)
entrypoint per pipeline:

| entrypoint | script | environment |
|---|---|---|
| `seurat` | `pipelines/seurat.R` | `envs/r` |
| `osca` | `pipelines/osca.R` | `envs/r` |
| `scrapper` | `pipelines/scrapper.R` | `envs/r` |
| `bixverse` | `pipelines/bixverse.R` | `envs/bixverse` |

Shared code lives in `src/`: argument parsing (`cli.R`), input loading and
output writing (`common.R`), the resolution search (`search.R`) and phase
timings (`timing.R`), which are recorded as [obkit](https://github.com/omnibenchmark/obkit)
events in `obkit-events.jsonl`.

## Environments

Each `envs/<env>/` is its own [rv](https://github.com/A2-ai/rv) project (R 4.5,
Bioconductor 3.22). CI (`.github/workflows/containers.yml`) builds each into
`ghcr.io/scrna-bench/pipelines-r:<env>-<sha>`, which the benchmark plan
references. bixverse builds from Rust sources, so its image is built with a
Rust toolchain that stays out of the final image.
