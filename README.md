# Santa's Stolen Sleigh

Code for Santa's Stolen Sleigh. Ranked 7th / xxxxx teams.

- Hill climbing method using *destroy and repair* operation.
- Fast *weighted trip length* function with insert operation. It can be computed in constant time.
- Many duplicate code and magic numbers. Sorry about that.

## Requirements

- g++ 4.7 or later
- make

# Compiling & Running

Place the [data files](https://www.kaggle.com/c/santas-stolen-sleigh/data) into a subfolder `./data`.

```
% ls data
gifts.csv
```

Compile.
```
% make
```

If you need specific version of g++,
```
% make clean
% GXX=g++-5.2 make
```

Run.
```
% ./santa
```

|Processing Time|Total weighted reindeer weariness|
|---------------|-------|
|2 hours| 12.4185B|
|4 hours| 12.4113B|
|8 hours| 12.4064B|
|16 hours| 12.4046B|
|32 hours| 12.4034B|
|42 hours| 12.4030B|
|...| ... |
|4.6 billion years| ? |

For more details, See [log.txt](log.txt).
