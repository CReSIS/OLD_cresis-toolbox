# Viterbi Costs

Documentation of the extent to which various factors influence binary and unary cost calculation used in finding the best path in [viterbi.cpp](viterbi.cpp).

## Context

The *unary cost* is calculated for a given (x, y) point without regard for surrounding points (except in the case of image magnitude correlation cost). This occurs within `double* viterbi::find_path(void)`.

The *binary cost* is calculated within [viterbi.h](viterbi.h) inside `void dt(...)` based on the distance from a given point to another.

The *cummulative cost* is stored in the arrays `path_prob` and `path_prob_next` and accumulates the unary and binary cost to each row for the current and previous columns.

Finding the best path is done by determining which row in the previous column to a given row in the current column yields the lowest cummulative cost. This is done for every row in the current column thus producing an array of lowest cost paths to every row in the current column. This is repeated for the next column and so on, only having to consider the previous column in each step (as the rows of the previous column are already known to hold the lowest cummulative cost possible to that point).

## Trace

Simplified steps taken in the Viterbi implementation.

- `viterbi(...)`  
  - calls `find_path()`  
    - initializes `path`, `path_prob`, `path_prob_next`, and `index` to 0 arrays of size `depth`  
    - calls `viterbi_right(...)`  
      - for every column:  
        - appends `index` array to `path` matrix if after first column
        - adds `unary_cost(col, row)` to every row
        - calculates scale (`norm`) if not last column    
        - calls `dt_1d(...)` (distance transform -- binary cost calculation) if not last column  
          - calls `dt(...)` with initial recursion conditions  
            - recursively populates `path_prob_next` with minimum cummulative cost to each row from the previous  
            - recursively populates `index` with previous row index corresponding to these costs  
        - swaps `path_prob` with `path_prob_next` to use current column costs as previous for next column
    - calls `calculate_best(path_prob)` and stores in `viterbi_index`
      - returns index of minimum cummulative cost in last column
    - populates `f_result` array back-to-front by following path of best indices backward from `viterbi_index` to the first column in `path`
    - returns `f_result`

## Cost Factors

Factor | Influence | Conditions | Implmentation | Notes | Upper end | Lower end
---|---|---|---|---|---|---
No Ice | `LARGE` | | `return LARGE`;
Above Surface | `LARGE` | | `return LARGE`;
Far from Center GT | `LARGE` | | `return LARGE`;
Far from Extra GT | Quadratic |
Surface Multiple Bin | Exponential |
Far from Ice Mask | Linear |
Image Magnitude | Negative Linear | | | Decreases Cost
Binary Cost (dt) | Quadratic |