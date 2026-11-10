| Fault | Faults | Expected | Detected | As expected | Caught by |
| :-- | :-- | :-- | --: | --: | :-- |
| Features deleted, duplicated or retyped in the database | F01, F02, F03, F09, F11, F19 | detected | 18 of 18 | 18 of 18 | source count, parent links |
| Parent links removed from the database | F04, F05, F18 | detected | 9 of 9 | 9 of 9 | parent links |
| Stored assembly damaged or moved | F13, F14, F15, F17 | detected | 12 of 12 | 12 of 12 | stored assembly |
| Input file edited after the load | F16 | detected | 3 of 3 | 3 of 3 | input's SHA-256 (count skipped, exit 0) |
| Parser drops, invents or retypes rows | P01, P02, P03, P04, P06 | detected | 15 of 15 | 15 of 15 | load refused: parent not defined, load rolled back: counts short, load rolled back: links short, load rolled back: row stored twice |
| Values changed, counts kept (database) | F06, F07, F08, F10, F12 | not detected | 0 of 15 | 15 of 15 | – |
| Values changed, counts kept (parser) | P05, P07, P08 | not detected | 0 of 9 | 9 of 9 | – |
