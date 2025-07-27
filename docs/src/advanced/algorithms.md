# [Algorithm Documentation](@id algorithm-documentation)

## Spline Evaluation

Fragment intensities are modeled with cubic B-splines.  Pioneer now
uses [De Boor's algorithm](https://en.wikipedia.org/wiki/De_Boor%27s_algorithm)
for evaluating these splines efficiently.  The implementation is found
in `src/utils/ML/libraryBSpline.jl` as `splevl_fast` and replaces the
older recursive evaluator.