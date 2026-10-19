# Changelog

## 0.1.0 (2026-10-19)


### Features

* adapt workflow for use with AoU ([7dbbc3b](https://github.com/gymrek-lab/happler/commit/7dbbc3b2879b3b2cd65e0e2da23230058767b5ab))
* convert perform_test to JAX with JIT compilation for 200x speedup ([71362a7](https://github.com/gymrek-lab/happler/commit/71362a758b92405bb5e6138575fdb12efb63592b))
* finish running `platelet_count` pipeline on UKB ([2bad3c3](https://github.com/gymrek-lab/happler/commit/2bad3c33f6b5d0b6c906ea3343dac234e8aa471e))
* iterative approach ([3090781](https://github.com/gymrek-lab/happler/commit/309078153c47db5b7ae2b31af7c09281f6acf3f4))
* MAF analysis (Fig 5) in AoU ([7668902](https://github.com/gymrek-lab/happler/commit/7668902bd5ed77b3f09c6b7f76e4efb8f85ce353))


### Bug Fixes

* apply new thresholds ([43d39d5](https://github.com/gymrek-lab/happler/commit/43d39d5f2eb4bb35f22058b57c4ffeaefe2c91f0))
* calibrate default thresholds ([91f25ab](https://github.com/gymrek-lab/happler/commit/91f25abe7a0f48777a785a617453d39ebd505d69))
* log transformations with bic in midway_manhattan_summary.py ([a682bc6](https://github.com/gymrek-lab/happler/commit/a682bc69179191f73d816523b3e9e3ed1a3da621))
* pgenlibr issue ([a506f08](https://github.com/gymrek-lab/happler/commit/a506f089f8eae4e8b6d55034fc36fdb842109e99))
* rename project from 'happler' to 'happler-fm' to resolve name collision on PyPI ([962921c](https://github.com/gymrek-lab/happler/commit/962921c011209c183ea45f31cf92e3273c66ffc1))


### Performance Improvements

* add original code for non-JAX FastBIC ([692640f](https://github.com/gymrek-lab/happler/commit/692640fb17e9008d48e351d5606d972ab63c4f8e))
* switch to ranking and thresholding only on BIC ([a1db690](https://github.com/gymrek-lab/happler/commit/a1db69066d1077489f3e106d62e680716c669d62))
* use JAX for 2x speedup with transform and sum operations ([deeda9b](https://github.com/gymrek-lab/happler/commit/deeda9b9ff842d9087bc40fb08d7f14ad2750774))


### Documentation

* describe happler commands ([53aa402](https://github.com/gymrek-lab/happler/commit/53aa402f68a0beec4fd8486f674fce5660b3e142))
