# ppm_pngc
PEP Mutase and PngC - code for Monika Papkinski's masters.

## TODO (specifically wrt. `protein_distances.py`)
- actually factor code into functions (parse_lastcol is some kind of a start)
- medium: fix exceptions
- medium: get it to only show the smallest distance per pep_mutase
- make it so that we don't need to refer to any columns by number (stretch goal :P)
- prefer/do pre-processing/filtering over post-processing - it's better and faster for us to compute less data overall (doi, i guess) even though it is annoying to preprocess some of these things...(or impossible, maybe)
- easy probably?: use tabs instead of generic spaces in output (??)
- catch more specific exceptions/get code to be safe enough to not need the big try/catch
- easy: turn the path into an optional command line argument. (that way, you can specify scaffold/complete on the command line instead of writing a new file for it *staring intently*)


## Dev & Debugging Notes/Tips
- use `ipdb` liberally! `ipython` is also a great tool.
- when testing, always use a subset of the total dataset. this can be as simple as limiting outfiles in that for loop like:
```
for file in out_files[:5]:
# instead of
for file in out_files:
```
- the folder is called virtual env, but we don't use virtual envs, really, so don't bother worrying about it. we do however have a requirements.txt which is decent enough. (I think I used this stuff because we had to install pip packages locally on the wlu servers)
- pandas has a bit of hard to learn syntax things, but my sentiment in using it was, "Oh, I can access columns by *name*". Unfortunately, none of the data files have properly formatted headers so I manually added one, perhaps adding the column names for GFF files is in our future...
- pandas has a LOT of utility functions. for example, `dataframe.drop_duplicates('column_name')` . their website has pretty good documentation.
- doing dev on the servers sucks (imo, so i never do it). However if you're not on the servers, note that you can't work remotely using SSH. => You *need* to have a test data set on your computer if you're doing dev on your own computer.
- please use something like pastebin or gist to post outputs to me. I *hate* looking at screenshots of code/output (unless it's funny)
- feel free to message me if you need help/have questions! this is far from perfect, or even good.
