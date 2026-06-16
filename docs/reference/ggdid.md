# Plot `did` objects using `ggplot2`

Function to plot objects from the `did` package

## Usage

``` r
ggdid(object, ...)
```

## Arguments

- object:

  either a `MP` object or `AGGTEobj` object. See
  [`help(ggdid.MP)`](https://bcallaway11.github.io/did/reference/ggdid.MP.md)
  and
  [`help(ggdid.AGGTEobj)`](https://bcallaway11.github.io/did/reference/ggdid.AGGTEobj.md).

- ...:

  other arguments

## Examples

``` r
if (FALSE) { # \dontrun{
data(mpdta)
out <- att_gt(yname = "lemp",
              gname = "first.treat",
              idname = "countyreal",
              tname = "year",
              xformla = ~1,
              data = mpdta)

# plot all group-time average treatment effects
ggdid(out)

# plot event study aggregation
es <- aggte(out, type = "dynamic")
ggdid(es)
} # }
```
