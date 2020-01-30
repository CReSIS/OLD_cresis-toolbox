# Changelog

In no particular order, every change made by me regarding Viterbi:

```markdown
Moved unary cost addition to before binary cost addition to account for first column unary cost
- Removed 2 column padding in main viterbi loop

Removed influence of mu in image magnitude correlation

Fixed EGT weight calculation (was always negative)
- Rounded EGT locations to ints to correspond with indexes of viterbi data

Set mu_size to 0 in picktool

Added scaling to picker viterbi_data to match echogram size in picktool
- Added scaling back to match echo size in picktool
- Interped gt to match y axis in picktool

Added multiple suppression to unary cost based on plane and surface position
- Added plane bin argument to viterbi constructor
- Added constants to viterbi.h to influence this effect on cost

Added viterbi_costs.md to document influence of various factors on cost
- Created this changelog

Corrected various error messages regarding viterbi input sizes

Added several todos to code

Removed top and multiple suppression from picktool
```
