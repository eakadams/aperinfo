# aperinfo
Collect information for Apertif imaging survey

Packages required:
- atdbquery (Note - updated to not return all status types to speed things up)

To update atdb obs list:<br>
`import info.atdb as atdb`<br>
`atdb.get_obstable()`<br>
Note - it's useful to do this in a screen as it takes quite a while to run.

