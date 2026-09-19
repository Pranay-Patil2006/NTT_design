# NTT Design

Hardware implementation of Number Theoretic Transform (NTT) using Verilog. This repository contains various architectural variants for computing NTT-based negacyclic convolution.

## Overview

The designs implement negacyclic convolution over GF(257). It computes `c = INTT( NTT(a) ⊙ NTT(b) )` where `⊙` is point-wise Montgomery multiplication. 
The transform uses Cooley-Tukey for the forward layers and Gentleman-Sande for the inverse layers.

### Key Parameters
- **Field:** GF(257)
- **Transform Size:** 32-point (32 buses)
- **Bus Width:** 9 bits

## Architecture Variants

The `src/` directory includes three distinct implementations targeting different trade-offs in area, performance, and latency:

1. **Optimized Combinational (`ntt_optimized_combinational.v`)**
   - Fully combinational NTT-based convolution.
   - High throughput but potentially longer critical path depending on synthesis tools.

2. **Single Delay Feedback Pipelined (`ntt_single_delay_feedback_pipelined.v`)**
   - Implements a pipelined architecture using Single Delay Feedback (SDF).
   - Suitable for streaming data and achieving higher clock frequencies.

3. **Time-Multiplexed 3-Cycle (`ntt_time_multiplexed_3_cycle.v`)**
   - A 3-phase time-multiplexed NTT convolution.
   - Computes `c = INTT( NTT(a) ⊙ NTT(b) )` using a single shared unified transform core over 3 clock cycles.
   - **Phase 0:** `NTT(a)`
   - **Phase 1:** `NTT(b)`
   - **Phase 2:** `INTT(dot_product)`
   - More area-efficient by reusing the transform core.

## Directory Structure

```text
.
├── docs/                 # Documentation and reports
│   └── NTT.pdf           # Detailed project report
├── src/                  # Verilog RTL source files
│   ├── ntt_optimized_combinational.v
│   ├── ntt_single_delay_feedback_pipelined.v
│   └── ntt_time_multiplexed_3_cycle.v
└── README.md
```

## Documentation

For a detailed explanation of the mathematics, architecture, and design decisions, please refer to the report located at [`docs/NTT.pdf`](docs/NTT.pdf).

## Author

- **Pranay Arvind Patil**
- Indian Institute of Technology, Gandhinagar
