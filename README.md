# OPNA_1 - Continued Fraction Generator
`./opna_1 pi --table --iterations 5 --precision 10`
```
╔═══════════╦════════════════╦═══════════════╦════════════════════╦════════════════════════╗
║ Iteration ║ Indices        ║ Fraction      ║ Evaluated fraction ║ Difference             ║
╠═══════════╬════════════════╬═══════════════╬════════════════════╬════════════════════════╣
║ 1         ║ [3]            ║ 3/1*          ║ 3.0000000000       ║  1.4159265359 * 10^-01 ║
║ 2         ║ [3,7]          ║ 22/7*         ║ 3.1428571429       ║ -1.2644892673 * 10^-03 ║
║ 3         ║ [3,7,15]       ║ 333/106*      ║ 3.1415094340       ║  8.3219627529 * 10^-05 ║
║ 4         ║ [3,7,15,1]     ║ 355/113*      ║ 3.1415929204       ║ -2.6676418906 * 10^-07 ║
║ 5         ║ [3,7,15,1,292] ║ 103993/33102* ║ 3.1415926530       ║  5.7789063439 * 10^-10 ║
╚═══════════╩════════════════╩═══════════════╩════════════════════╩════════════════════════╝
```
```
Usage: ./opna_1 [OPTIONS] <NUMBER|CONSTANT> 
Allowed options:
  --help                 print help message
  --version              print version information
  --examples             show examples
  --table                print evaluation table of every iteration
  --table_style arg      Specify table style: nice|double|simple|empty
  --headerless           Do not display header when showing results
  --fields arg           Comma delimited list of fields to be shown: 
                         iter,ind,frac,eval,diff
  --maxdenominator arg   maximum denominator up to which to iterate
  --inbetween            show best in-between continual fraction approximations
                         as well
  --iterations arg (=14) number of iterations
  --precision arg (=14)  how many decimals should result have
  --number arg           number to be parsed

Allowed constants: pi, phi, e 
```
### Examples
```
    Command                                                Description                                                                                                  
 1  ./opna_1 pi                                            Evaluates pi with default (14) iterations                                                                    
 2  ./opna_1 pi --table --iterations 5 --precision 10      Print evaluation table with 5 iterations and up to 10 decimal places for pi expansion                        
 3  ./opna_1 phi --table --maxdenominator 400 --inbetween  Print evaluation table for phi where maximum fraction approximation denominator is less than or equal to 400 
 4  ./opna_1 phi --table --inbetween                       Print evaluation table for phi with default (14) iterations and also find best in-between approximations     
```

### Build requirements
 C++17 compatible compiler, boost libraries, cmake
