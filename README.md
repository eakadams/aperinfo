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


## DR2 notes
### Polarization
- Collect all new polarization validation from `/home/adebahr/apercal/ipython-notebooks/commissioning/polarisation/Polarisation_QA/polarisation_all/` (See helpful rsync command in combine_pol doc)
- Run `valid.combine_pol()`
- Run `valid.do_pol_valid()` to ensure consistency with before
  - Note that I need to check `do_pol_valid` for how taskids were manually set. This could be an appropriate place to manually set taskids that I know have bad pol data products

### Line
- Had to get cube status explicitly: wrote valid.update_line_valid for this:
  ```from info import valid
  valid.update_line_valid()```