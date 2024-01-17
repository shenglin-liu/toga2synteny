# toga2synteny
Generate a chromosomal synteny map among a group of species based on TOGA annotation

## Usage
Can be found at the beginning of the file.

## To do
Automate the plotting.

## Example plot
<img width="1217" alt="image" src="https://github.com/shenglin-liu/toga2synteny/assets/9432577/ea74342c-9ef5-4402-bcfa-7ab89554fd3e">
The plot on the left is the default plotting mode. The one on the right rearranges the chromosomes (also flips chromosomes whenever necessary; see the grey ones) for a better visualization.

## Advantages
No reliance on external software or any packages.
Fast.
Wrapped around TOGA results; but the core functions have more generic appliances.
Figures have small file sizes.

## Requirements/limitations
Good assembly quality, preferably chromosome level (N50 at least 10Mb).
Good annotation quality, not too many missing (only one2one orthologous genes existing in all species are used).
