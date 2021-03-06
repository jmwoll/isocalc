With isocalc, isotopic patterns can be calculated
from python. It is not accurate (e.g. not considering electrons),
and very inefficient, but still sometimes useful for quickly screening
isotopic patterns and for teaching. For example, the code
```
from isocalc import isocalc
isocalc.plot_isotope_distribution("C1")
```
will result in the following plot:

![Example Isotopic Pattern.](https://github.com/jmwoll/isocalc/blob/master/doc/isocalc_example_c.png)

While the following code:
```
from isocalc import isocalc
isocalc.plot_isotope_distribution("C1Br4")
```
will result in the following plot with several isotopic peaks:

![Example Isotopic Pattern.](https://github.com/jmwoll/isocalc/blob/master/doc/isocalc_example_cbr4.png)
