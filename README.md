# Global Basin

## LOWRES ##

I discovered a disturbing sensitivity of the solution to viscosity. In the original config, there was a huge amount of cooling and deep convection in the Southern Ocean. I noticed lots of grid-scale noise in W and decided to turn the viscosity from viscAhGrid=2.E-2 to viscAhGrid=1.E-1. The sign of the heat flux completely changed sign! There is still grid-scale noise, so maybe it has to be decreased further.
