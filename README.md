# PowerCompareSC

Semi-continuous data is characterized by a population that mostly normal/continuous except for a significant subset of 0 values. A an example would be number of cigarettes smoked per week by Americans. The data for smokers would clearly be continuos and likely normally distributed, however there would a significant percent of 0 values in the dataset representing non-smokers. This type of data, also often called zero-inflated data, is notoriously hard to properly analyze. Depending on the desired output and/or characterization of the specific dataset itself, one might use a non-parametric test such as the Wilcoxon test, or simply accept the caveats and choose to rely on a standard t-test.

This R package presents a new methodology for analysis of such zero-inflated data presenting an analysis method and test that will often times prove be more reliable and to provide more power under the same conditions. The package itself includes various functions to compare, test, and run new inference procedures or analysis of semi-continuous/zero-inflated data. 
 
