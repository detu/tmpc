#!/bin/bash
build/bin/tmpc_bench --benchmark_filter=.+Riccati --benchmark_out=benchmark.json --benchmark_repetitions=10 --benchmark_report_aggregates_only=true
