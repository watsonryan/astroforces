# Drag Performance Baseline

## Build
```bash
cmake --preset macos-debug
cmake --build --preset macos-debug -j
```

## Run
```bash
./build/macos-debug/astroforces_perf_benchmark --samples 40 --iters 5000
```

Optional controls:
```bash
./build/macos-debug/astroforces_perf_benchmark --samples 80 --iters 10000
```

## Metrics
- `single_eval_mean_us`: mean microseconds per single drag evaluate call
- `single_eval_hz`: equivalent single-call throughput
- `batch_eval_mean_us`: mean microseconds per sample in batch loop
- `batch_eval_hz`: equivalent per-sample throughput in batch mode
