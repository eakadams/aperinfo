# aperinfo
Collect information for Apertif imaging survey

Packages required:
- atdbquery (Note - updated to not return all status types to speed things up)
- kapteyn (For skyplots)

To update atdb obs list:<br>
`import info.atdb as atdb`<br>
`atdb.get_obstable()`<br>
Note - it's useful to do this in a screen as it takes quite a while to run.

To get a summary of observations to date:<br>
```
import overview as ov
obs = ov.ObsCat()
#print a summary of fields observed
obs.get_summary_obs()
#make a sky plot
obs.plot_all_obs()
```


